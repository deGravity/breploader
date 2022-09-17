#ifndef BODY_H_INCLUDED
#define BODY_H_INCLUDED 1

#include <parasolid.h>
#include <TopoDS_Shape.hxx>
#include "face.h"
#include "loop.h"
#include "edge.h"
#include "topology.h"


#include <map>
#include <unordered_map>
#include <Eigen/Core>
#include <iostream>
#include <memory>
#include "types.h"

namespace pspy {

class Body {
public:
    virtual BREPTopology GetTopology() = 0;

    virtual MassProperties GetMassProperties(double accuracy = MASS_ACC) = 0;

    virtual Eigen::MatrixXd GetBoundingBox() = 0;

    virtual int Transform(const Eigen::MatrixXd& xfrm) = 0;

    virtual void Tesselate(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::VectorXi& FtoT,
        Eigen::MatrixXi& EtoT,
        Eigen::VectorXi& VtoT) = 0;

    virtual void debug() = 0;
};

class PSBody: public Body {
public:
    PSBody(int id);
    ~PSBody();

    BREPTopology GetTopology();

    MassProperties GetMassProperties(double accuracy = MASS_ACC);

    Eigen::MatrixXd GetBoundingBox();

    int Transform(const Eigen::MatrixXd& xfrm);

    void Tesselate(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::VectorXi& FtoT,
        Eigen::MatrixXi& EtoT,
        Eigen::VectorXi& VtoT);

    void debug();

private:
    int _id;
    bool _valid; // If we need to re-compute due to transforms
};

class OCCTBody: public Body {
public:
    OCCTBody(const TopoDS_Shape& shape);

    BREPTopology GetTopology();

    MassProperties GetMassProperties(double accuracy = MASS_ACC);

    Eigen::MatrixXd GetBoundingBox();

    int Transform(const Eigen::MatrixXd& xfrm);

    void Tesselate(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::VectorXi& FtoT,
        Eigen::MatrixXi& EtoT,
        Eigen::VectorXi& VtoT);

    void debug();

private:
    TopoDS_Shape _shape;
    std::unordered_map<
        TopoDS_Shape,
        int,
        TopoDS_Shape_Hash<TopoDS_Shape>,
        TopoDS_Shape_Pred<TopoDS_Shape> > _shape_to_idx;
    bool _valid; // If we need to re-compute due to transforms
};

// Helper Functions
std::vector<std::shared_ptr<Body> > read_file(std::string path);

// PSBody Helper Functions
bool is_psbody(int id);
std::vector<std::shared_ptr<Body> > read_xt(std::string path);

// OCCTBody Helper Functions
std::vector<std::shared_ptr<Body> > read_step(std::string path);

}

#endif // !BODY_H_INCLUDED
