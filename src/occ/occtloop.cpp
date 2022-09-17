#include <loop.h>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepTools.hxx>
#include <ShapeExtend_WireData.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <Bnd_OBB.hxx>
#include <BRepBndLib.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <gp_Pln.hxx>

namespace pspy {

OCCTLoop::OCCTLoop(const TopoDS_Shape& shape, const TopTools_ListOfShape& faces)
{
    assert(shape.ShapeType() == TopAbs_WIRE);
    
    _shape = TopoDS::Wire(shape);
    _faces = faces;

    // Set loop to be inner if it is not the outer loop of one of
    // its faces
    _type = LoopType::OUTER;

    for (auto iterator = faces.cbegin(); iterator != faces.cend(); iterator++) {
        TopoDS_Face face_shape = TopoDS::Face(*iterator);
        TopoDS_Wire face_outer_wire = BRepTools::OuterWire(face_shape);

        if (face_outer_wire.IsNotEqual(_shape)) {
            _type = LoopType::INNER;
        }
    }

    ShapeExtend_WireData loop_data = ShapeExtend_WireData(_shape);
    int num_edges = loop_data.NbEdges();

    // A loop can be degenerate (e.g. the tip of a cone)
    // so we can't always get mass properties
    if (num_edges > 0) {
        // Check if this loop is just a single circle
        // this is needed in inference calculation
        // TODO: fix for loops with multiple edges that are still just a single circle
        if (num_edges == 1) {
            TopoDS_Edge edge = loop_data.Edge(1);
            BRepAdaptor_Curve curve(edge);

            if (curve.GetType() == GeomAbs_Circle) {
                _is_circle = true;
            }
            else {
                _is_circle = false;
            }
        }
        else {
            _is_circle = false;
        }

        auto m = MassProperties(_shape);
        length = m.amount;
        center_of_gravity = m.c_of_g;
        moment_of_inertia = m.m_of_i;

        // Get Non-Aligned Bounding Box
        // Axes are ordered so X is largest, Y is second largest, Z is smallest
        // TODO: check for null curves?
        Bnd_OBB nabox;
        BRepBndLib::AddOBB(_shape, nabox);

        gp_XYZ center = nabox.Center();
        gp_XYZ xdirection = nabox.XDirection();
        gp_XYZ zdirection = nabox.ZDirection();
        double xhsize = nabox.XHSize();
        double yhsize = nabox.YHSize();
        double zhsize = nabox.ZHSize();

        na_bb_center <<
            center.X(),
            center.Y(),
            center.Z();
        na_bb_x <<
            xdirection.X(),
            xdirection.Y(),
            xdirection.Z();
        na_bb_z <<
            zdirection.X(),
            zdirection.Y(),
            zdirection.Z();
        na_bounding_box.resize(2, 3);
        na_bounding_box <<
            -xhsize, -yhsize, -zhsize,
            xhsize, yhsize, zhsize;
    }
    else {
        length = 0;
        center_of_gravity = Eigen::Vector3d(0, 0, 0);
        moment_of_inertia = Eigen::MatrixXd::Zero(3, 3);

        // We could techincally get this for vertices if they exist
        na_bb_center = Eigen::Vector3d(0, 0, 0);
        na_bb_x = Eigen::Vector3d(0, 0, 0);
        na_bb_z = Eigen::Vector3d(0, 0, 0);
        na_bounding_box = Eigen::MatrixXd::Zero(2, 3);
    }
}

std::vector<Inference> OCCTLoop::get_inferences()
{
    std::vector<Inference> inferences;

    // Only make inferences for inner loops
    if (_type == LoopType::INNER) {
        // Make inferences for plane for which loop is an inner loop
        for (auto iterator = _faces.cbegin(); iterator != _faces.cend(); iterator++) {
            TopoDS_Face face_shape = TopoDS::Face(*iterator);
            BRepAdaptor_Surface face_surf(face_shape);
            TopoDS_Wire face_outer_wire = BRepTools::OuterWire(face_shape);

            if (face_surf.GetType() == GeomAbs_Plane && face_outer_wire.IsNotEqual(_shape)) {
                gp_Pln plane = face_surf.Plane();
                TopAbs_Orientation face_orientation = face_shape.Orientation();

                // Get the plane's axis to use
                Eigen::Vector3d plane_axis(
                    plane.Axis().Direction().X(),
                    plane.Axis().Direction().Y(),
                    plane.Axis().Direction().Z());
                if (face_orientation == TopAbs_REVERSED) {
                    plane_axis = -plane_axis;
                }

                Inference inf;
                inf.inference_type = InferenceType::LOOP_CENTER;
                inf.origin = center_of_gravity;
                inf.z_axis = plane_axis;

                // Onshape Special Case - Onshape does not allow
                // loop center inferences for circular loops
                if (_is_circle) {
                    inf.onshape_inference = false;
                }

                inferences.push_back(inf);
            }
        }
    }

    return inferences;
}

}