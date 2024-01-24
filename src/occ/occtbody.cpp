#include <body.h>

#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>

#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopExp.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>
#include <Interface_Static.hxx>
#include <gp_Trsf.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopExp_Explorer.hxx>

namespace pspy {

std::vector<std::shared_ptr<Body>> read_step(std::string path) {
    std::vector<std::shared_ptr<Body>> parts_vec;

	STEPControl_Reader reader;
    Interface_Static::SetCVal("xstep.cascade.unit", "M");
	IFSelect_ReturnStatus ret;

    // Allow the file contents to be passed in instead of the path
    // All step files begin with "ISO-10303-21;", and this is not
    // a valid file path, so this should be safe.
    if (path.rfind("ISO-10303-21;", 0) == std::string::npos) {
        ret = reader.ReadFile(path.c_str());
    } else {
        auto path_ss = std::stringstream(path);
        ret = reader.ReadStream("rawtext", path_ss);
    }
    if (ret == IFSelect_RetDone) {
        auto num_roots = reader.TransferRoots();
        if (num_roots == 1) {
            TopoDS_Shape shape = reader.OneShape();
            parts_vec.emplace_back(new OCCTBody(shape));
        }
    }

    return parts_vec;
}

OCCTBody::OCCTBody(const TopoDS_Shape& shape) {
    _shape = shape;

    // Create _shape_to_idx to associate entities with logical IDs
    TopTools_IndexedMapOfShape shape_map;
    TopExp::MapShapes(_shape, shape_map);

    int i = 1;
    for (auto iterator = shape_map.cbegin(); iterator != shape_map.cend(); iterator++) {
        TopoDS_Shape subshape = *iterator;
        _shape_to_idx[subshape] = i;
        i++;
    }

    _valid = true;
}

BREPTopology OCCTBody::GetTopology() {
    BREPTopology topology;

    std::map<int, int> cat_idx;

    TopTools_IndexedDataMapOfShapeListOfShape edge_face_map;
    TopExp::MapShapesAndAncestors(_shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map);

    TopTools_IndexedDataMapOfShapeListOfShape loop_face_map;
    TopExp::MapShapesAndAncestors(_shape, TopAbs_WIRE, TopAbs_FACE, loop_face_map);

    // Store faces, loops, edges, vertices
    // Iterate in order to maintain ID stability
    TopTools_IndexedMapOfShape shape_map;
    TopExp::MapShapes(_shape, shape_map);

    int i = 0;
    for (auto iterator = shape_map.cbegin(); iterator != shape_map.cend(); iterator++) {
        TopoDS_Shape subshape = *iterator;
        int idx = _shape_to_idx[subshape];

        switch (subshape.ShapeType()) {
        case TopAbs_FACE:
            topology.pk_to_class[idx] = TopologyType::FACE;
            cat_idx[i] = topology.faces.size();
            topology.faces.emplace_back(new OCCTFace(subshape));
            break;
        case TopAbs_WIRE:
            topology.pk_to_class[idx] = TopologyType::LOOP;
            cat_idx[i] = topology.loops.size();
            topology.loops.emplace_back(new OCCTLoop(subshape, loop_face_map.FindFromKey(subshape)));
            break;
        case TopAbs_EDGE:
            topology.pk_to_class[idx] = TopologyType::EDGE;
            cat_idx[i] = topology.edges.size();
            topology.edges.emplace_back(new OCCTEdge(subshape, edge_face_map.FindFromKey(subshape)));
            break;
        case TopAbs_VERTEX:
            topology.pk_to_class[idx] = TopologyType::VERTEX;
            cat_idx[i] = topology.vertices.size();
            topology.vertices.emplace_back(new OCCTVertex(subshape));
            break;
        default:
            break;
        }

        // Update map from pk index to index with each entity type
        topology.pk_to_idx[idx] = cat_idx[i];

        i++;
    }

    // Map for finding face-face edges
    std::map<int, std::vector<int> > edge_to_faces;

    std::map<int, std::vector<int> > face_to_loops;
    std::map<int, std::vector<int> > face_to_edges;
    std::map<int, std::vector<int> > face_to_vertices;
    std::map<int, std::vector<int> > loop_to_edges;
    std::map<int, std::vector<int> > loop_to_vertices;
    std::map<int, std::vector<int> > edge_to_vertices;

    // Fill topology maps
    for (const auto& subshape_idx_pair : _shape_to_idx) {
        TopoDS_Shape subshape = subshape_idx_pair.first;
        int idx = subshape_idx_pair.second;
        int parent = topology.pk_to_idx[idx];
        int child;

        if (subshape.ShapeType() == TopAbs_FACE) {
            // Fill wires
            TopExp_Explorer explorer(subshape, TopAbs_WIRE);
            while (explorer.More()) {
                TopoDS_Shape wire = explorer.Current();
                child = topology.pk_to_idx[_shape_to_idx[wire]];

                topology.face_to_loop.emplace_back(parent, child, TopoRelationSense::None);
                face_to_loops[parent].push_back(child);

                explorer.Next();
            }

            // Fill edges
            explorer.Init(subshape, TopAbs_EDGE);
            while (explorer.More()) {
                TopoDS_Shape edge = explorer.Current();
                child = topology.pk_to_idx[_shape_to_idx[edge]];

                face_to_edges[parent].push_back(child);
                edge_to_faces[child].push_back(parent);

                explorer.Next();
            }

            // Fill vertices
            explorer.Init(subshape, TopAbs_VERTEX);
            while (explorer.More()) {
                TopoDS_Shape vertex = explorer.Current();
                child = topology.pk_to_idx[_shape_to_idx[vertex]];

                face_to_vertices[parent].push_back(child);

                explorer.Next();
            }

        }
        else if (subshape.ShapeType() == TopAbs_WIRE) {
            // Fill edges
            BRepTools_WireExplorer wire_explorer(TopoDS::Wire(subshape));
            while (wire_explorer.More()) {
                TopoDS_Edge edge = wire_explorer.Current();
                TopAbs_Orientation orientation = wire_explorer.Orientation();
                child = topology.pk_to_idx[_shape_to_idx[edge]];
                TopoRelationSense sense = 
                    orientation == TopAbs_FORWARD ? TopoRelationSense::Positive :
                    orientation == TopAbs_REVERSED ? TopoRelationSense::Negative :
                    TopoRelationSense::None;

                topology.loop_to_edge.emplace_back(parent, child, sense);
                loop_to_edges[parent].push_back(child);

                wire_explorer.Next();
            }

            // Fill vertices
            TopExp_Explorer explorer(subshape, TopAbs_VERTEX);
            while (explorer.More()) {
                TopoDS_Shape vertex = explorer.Current();
                child = topology.pk_to_idx[_shape_to_idx[vertex]];

                loop_to_vertices[parent].push_back(child);

                explorer.Next();
            }
        }
        else if (subshape.ShapeType() == TopAbs_EDGE) {
            // Fill vertices
            TopExp_Explorer explorer(subshape, TopAbs_VERTEX);
            while (explorer.More()) {
                TopoDS_Shape vertex = explorer.Current();
                child = topology.pk_to_idx[_shape_to_idx[vertex]];

                topology.edge_to_vertex.emplace_back(parent, child, TopoRelationSense::None);
                edge_to_vertices[parent].push_back(child);

                explorer.Next();
            }
        }
    }

    // Also find Face-Face Edges
    for (auto edgefaces : edge_to_faces) {
        int edge = edgefaces.first;
        auto faces = edgefaces.second;
        if (faces.size() > 0) {
            int face1 = faces[0];
            int face2 = faces[faces.size() - 1];
            topology.face_to_face.emplace_back(face1, face2, edge);
        }
    }

    // Assign to structure
    // TODO - don't use the temporary variables
    // TODO - int -> size_t
    topology.face_loop = face_to_loops;
    topology.face_edge = face_to_edges;
    topology.face_vertex = face_to_vertices;
    topology.loop_edge = loop_to_edges;
    topology.loop_vertex = loop_to_vertices;
    topology.edge_vertex = edge_to_vertices;

    return topology;
}

MassProperties OCCTBody::GetMassProperties(double accuracy) {
    return MassProperties(_shape, accuracy);
}

Eigen::MatrixXd OCCTBody::GetBoundingBox() {
    Bnd_Box box;
    BRepBndLib::Add(_shape, box);

    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    Eigen::MatrixXd corners(2, 3);
    corners <<
        xmin, ymin, zmin,
        xmax, ymax, zmax;
    
    return corners;
}

int OCCTBody::Transform(const Eigen::MatrixXd& xfrm) {
    // TODO: perspective division, if desired.
    //  We currently assume the last row of xfrm is just
    //  (0, 0, 0, k) for some constant k.
    // Apply transform.
    double k = xfrm(3, 3);
    gp_Trsf transform_mat;
    transform_mat.SetValues(
        xfrm(0, 0) / k, xfrm(0, 1) / k, xfrm(0, 2) / k, xfrm(0, 3) / k,
        xfrm(1, 0) / k, xfrm(1, 1) / k, xfrm(1, 2) / k, xfrm(1, 3) / k,
        xfrm(2, 0) / k, xfrm(2, 1) / k, xfrm(2, 2) / k, xfrm(2, 3) / k
    );
    BRepBuilderAPI_Transform transform(_shape, transform_mat);

    // Remake _shape and _shape_to_idx.
    _shape = transform.ModifiedShape(_shape);
    auto shape_to_idx = _shape_to_idx;
    _shape_to_idx.clear();
    for (const auto& subshape_idx_pair : shape_to_idx) {
        TopoDS_Shape subshape = subshape_idx_pair.first;
        int idx = subshape_idx_pair.second;

        _shape_to_idx[transform.ModifiedShape(subshape)] = idx;
    }

    return 0;
}

void OCCTBody::Tesselate(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::VectorXi& FtoT,
    Eigen::MatrixXi& EtoT,
    Eigen::VectorXi& VtoT,
    bool set_quality,
    double quality) {
    // Setup faceting call options
    // TODO: give control over linear deflection
    double linear_deflection = 0.01;

    if (set_quality) {
        linear_deflection = quality;
    }

    // LinearDeflection - limits distance between a curve and its tesselation
    // AngularDeflection limits angle between subsequent segments in a polyline
    // DeflectionInterior - limits distance between triangles and face interior
    // AngleInterior (B-spline faces only) limits angle between triangle corner normals not on face boundaries
    // face boundaries use angular deflection instead

    // Facet the body
    BRepMesh_IncrementalMesh(_shape, linear_deflection ,false, 0.5, true);

    // Gather vertices and triangles from faces of _shape
    std::unordered_map<gp_Pnt, int, gp_Pnt_Hash, gp_Pnt_Pred> pnt_idxs;
    std::vector<gp_Pnt> pnts;
    std::vector<Eigen::Vector3i> tris;
    std::vector<int> tri_face;
    std::map<std::pair<int, int>, int> fin_edge;
    std::map<int, int> pnt_vertex;

    TopLoc_Location loc;

    // Fill pnt_idxs, pnts, tris, tri_face, and fin_edge
    for (const auto& subshape_idx_pair : _shape_to_idx) {
        TopoDS_Shape subshape = subshape_idx_pair.first;
        int subshape_idx = subshape_idx_pair.second;

        if (subshape.ShapeType() == TopAbs_FACE) {
            TopoDS_Face subface = TopoDS::Face(subshape);
            auto subface_triangulation =
                BRep_Tool::Triangulation(subface, loc);
            
            if (!subface_triangulation.IsNull()) {
                // Add new points
                for (int i = 1; i <= subface_triangulation->NbNodes(); ++i) {
                    gp_Pnt pnt = subface_triangulation->Node(i);
                    if (pnt_idxs.find(pnt) == pnt_idxs.end()) {
                        pnt_idxs[pnt] = pnts.size();
                        pnts.push_back(pnt);
                    }
                }

                // Add new triangles
                // Associate each fin (represented by a start point index and end point index
                // in the triangulation nodes) with a tri index in tris and the index
                // of the corresponding fin in the tri to allow filling
                // of fin_edge
                std::map<std::pair<int, int>, std::pair<int, int>> nodes_to_tri_fin;

                for (int i = 1; i <= subface_triangulation->NbTriangles(); ++i) {
                    Poly_Triangle tri = subface_triangulation->Triangle(i);
                    int node1 = tri.Value(1);
                    int node2 = tri.Value(2);
                    int node3 = tri.Value(3);
                    gp_Pnt pnt1 = subface_triangulation->Node(node1);
                    gp_Pnt pnt2 = subface_triangulation->Node(node2);
                    gp_Pnt pnt3 = subface_triangulation->Node(node3);
                    int pnt1_idx = pnt_idxs[pnt1];
                    int pnt2_idx = pnt_idxs[pnt2];
                    int pnt3_idx = pnt_idxs[pnt3];
                    int tri_idx = tris.size();
                    tris.emplace_back(pnt1_idx, pnt2_idx, pnt3_idx);
                    tri_face.push_back(subshape_idx);
                    nodes_to_tri_fin[std::pair<int, int>(node1, node2)] =
                        std::pair<int, int>(tri_idx, 1);
                    nodes_to_tri_fin[std::pair<int, int>(node2, node1)] =
                        std::pair<int, int>(tri_idx, 1);
                    nodes_to_tri_fin[std::pair<int, int>(node2, node3)] =
                        std::pair<int, int>(tri_idx, 2);
                    nodes_to_tri_fin[std::pair<int, int>(node3, node2)] =
                        std::pair<int, int>(tri_idx, 2);
                    nodes_to_tri_fin[std::pair<int, int>(node3, node1)] =
                        std::pair<int, int>(tri_idx, 0);
                    nodes_to_tri_fin[std::pair<int, int>(node1, node3)] =
                        std::pair<int, int>(tri_idx, 0);
                }

                // Add new fins associated with edges
                // Create _shape_to_idx to associate entities with logical IDs
                TopExp_Explorer explorer(subface, TopAbs_EDGE);
                while (explorer.More()) {
                    TopoDS_Edge edge = TopoDS::Edge(explorer.Current());
                    int edge_idx = _shape_to_idx[edge];

                    auto poly = BRep_Tool::PolygonOnTriangulation(edge, subface_triangulation, loc);
                    if (!poly.IsNull()) {
                        // Loop through the polygon, matching each fin with its index in fin_indices
                        // and matching that index with the index of the edge in _shape_to_idx
                        for (int i = poly->Nodes().Lower(); i < poly->Nodes().Upper(); ++i) {
                            int node1 = poly->Nodes()[i];
                            int node2 = poly->Nodes()[i + 1];

                            auto fin_idx = nodes_to_tri_fin[std::pair<int, int>(node1, node2)];
                            fin_edge[fin_idx] = edge_idx;
                        }
                    }

                    explorer.Next();
                }
            }
        }
    }

    // Fill pnt_vertex
    for (const auto& subshape_idx_pair : _shape_to_idx) {
        TopoDS_Shape subshape = subshape_idx_pair.first;
        int subshape_idx = subshape_idx_pair.second;

        if (subshape.ShapeType() == TopAbs_VERTEX) {
            TopoDS_Vertex subvertex = TopoDS::Vertex(subshape);
            gp_Pnt subvertex_pnt = BRep_Tool::Pnt(subvertex);

            if (pnt_idxs.find(subvertex_pnt) != pnt_idxs.end()) {
                pnt_vertex[pnt_idxs[subvertex_pnt]] = subshape_idx;
            }
        }
    }

    // Populate Mesh Vertices
    V.resize(pnts.size(), 3);
    for (int i = 0; i < pnts.size(); ++i) {
        auto pnt = pnts[i];
        V(i, 0) = pnt.X();
        V(i, 1) = pnt.Y();
        V(i, 2) = pnt.Z();
    }

    // Populate Mesh Faces
    F.resize(tris.size(), 3);
    for (int i = 0; i < tris.size(); ++i) {
        auto vec = tris[i];
        F(i, 0) = vec(0);
        F(i, 1) = vec(1);
        F(i, 2) = vec(2);
    }

    // Populate Mesh to Topology Face References
    FtoT.resize(tri_face.size());
    for (int i = 0; i < tri_face.size(); ++i) {
        FtoT[i] = tri_face[i];
    }
    // Populate Mesh to Topology Edge References
    EtoT = Eigen::MatrixXi::Zero(F.rows(), 3);
    for (const auto& fin_idx_edge_idx : fin_edge) {
        auto fin_idx = fin_idx_edge_idx.first;
        auto edge_idx = fin_idx_edge_idx.second;
        EtoT(fin_idx.first, fin_idx.second) = edge_idx;
    }
    // Populate Mesh to Topology Vertex References
    VtoT = Eigen::VectorXi::Zero(V.rows());
    for (const auto& pnt_idx_vertex_idx : pnt_vertex) {
        auto pnt_idx = pnt_idx_vertex_idx.first;
        auto vertex_idx = pnt_idx_vertex_idx.second;
        VtoT(pnt_idx) = vertex_idx;
    }
}

void OCCTBody::debug() {
    int err = 0; // PK_ERROR_t err = PK_ERROR_no_errors;
    TopAbs_ShapeEnum entity_class = _shape.ShapeType();
    std::cout << entity_class << std::endl;
    auto topo = GetTopology();
    auto mass = GetMassProperties();

    Eigen::MatrixXd xfrm(4, 4);
    xfrm <<
        1, 0, 0, 0.0,
        0, 1, 0, 0.0,
        0, 0, 1, 0.0,
        0, 0, 0, 0.5;
    err = Transform(xfrm);
    std::cout << "xfrm error = " << err << std::endl;
    auto topo_xfrmed = GetTopology();
    auto mass_xfrmed = GetMassProperties();


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi FtoT;
    Eigen::MatrixXi EtoT;
    Eigen::VectorXi VtoT;
    Tesselate(V, F, FtoT, EtoT, VtoT);

    std::cout << "Num Vertices: " << topo.vertices.size() << std::endl;
    std::cout << "Num face-face relations: " << topo.face_to_face.size() << std::endl;
    std::cout << "Num Tess Verts: " << V.rows() << std::endl;
    std::cout << "Num Tess Faces: " << F.rows() << std::endl;
    std::cout << "Mass Amount: " << mass.amount << std::endl;
    std::cout << "Mass Mass: " << mass.mass << std::endl;
    std::cout << "Mass c of G:" << std::endl << mass.c_of_g << std::endl;
    std::cout << "Mass m of i:" << std::endl << mass.m_of_i << std::endl;
    std::cout << "Num Tess Faces: " << F.rows() << std::endl;
}

}
