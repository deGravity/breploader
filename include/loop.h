#ifndef LOOP_H_INCLUDED
#define LOOP_H_INCLUDED 1

#include <parasolid.h>
#include "types.h"
#include <assert.h>
#include <Eigen/Core>
#include <vector>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListOfShape.hxx>

namespace pspy {

struct Loop {
    virtual std::vector<Inference> get_inferences() = 0;

    LoopType _type;

    bool _is_circle;

    double length;
    Eigen::Vector3d center_of_gravity;
    Eigen::MatrixXd moment_of_inertia;

    Eigen::Vector3d na_bb_center;
    Eigen::Vector3d na_bb_x;
    Eigen::Vector3d na_bb_z;
    Eigen::MatrixXd na_bounding_box;
};

struct PSLoop: public Loop {
    PSLoop(int id);

    std::vector<Inference> get_inferences() override;

    int _id;
};

struct OCCTLoop: public Loop {
    OCCTLoop(const TopoDS_Shape& shape, const TopTools_ListOfShape& faces);

    std::vector<Inference> get_inferences() override;

    TopoDS_Wire _shape;
    TopTools_ListOfShape _faces;
};

}

#endif // LOOP_H_INCLUDED
