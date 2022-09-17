#include <vertex.h>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>

namespace pspy {

OCCTVertex::OCCTVertex(const TopoDS_Shape& shape) {
    assert(shape.ShapeType() == TopAbs_VERTEX);

    _shape = TopoDS::Vertex(shape);

    gp_Pnt point = BRep_Tool::Pnt(_shape);

    position = Eigen::Vector3d(point.X(), point.Y(), point.Z());
}

std::vector<Inference> OCCTVertex::get_inferences()
{
    std::vector<Inference> inferences;

    Inference inf;
    inf.inference_type = InferenceType::POINT;
    inf.origin = position;
    inf.z_axis = Eigen::Vector3d(0.0, 0.0, 1.0);

    inferences.push_back(inf);

    return inferences;
}

}