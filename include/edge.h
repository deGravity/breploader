#ifndef EDGE_H_INCLUDED
#define EDGE_H_INCLUDED 1

#include "types.h"
#include <Eigen/Core>
#include <vector>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepAdaptor_Curve.hxx>

namespace pspy {

struct Edge {
    virtual std::vector<Inference> get_inferences() = 0;

    CurveFunction function;
    std::vector<double> parameters;

    double t_start;
    double t_end;

    Eigen::Vector3d start;
    Eigen::Vector3d end;
    bool _has_curve;
    bool _is_reversed; // Is edge opposite curve? Already baked-in to start and
                       // end, but needed for tangent computation.

    bool is_periodic;
    Eigen::Vector3d mid_point;

    double length;
    Eigen::Vector3d center_of_gravity;
    Eigen::MatrixXd moment_of_inertia;

    Eigen::MatrixXd bounding_box;
    Eigen::Vector3d na_bb_center;
    Eigen::Vector3d na_bb_x;
    Eigen::Vector3d na_bb_z;
    Eigen::MatrixXd na_bounding_box;

    virtual void sample_points(
        const int num_points,
        const bool sample_tangents,
        std::vector<Eigen::VectorXd>& samples,
        Eigen::Vector2d& t_range) = 0;

    virtual bool sample_curve(
        const int N_samples, // N
        Eigen::Vector2d& t_bounds,
        Eigen::MatrixXd& t_samples // (Nx6) x,y,z,t_x,t_y,t_z
    ) = 0;
};

struct PSEdge: public Edge {
    PSEdge(int id);

    void init_bb();
    void init_nabb();

    void init_line();
    void init_circle();
    void init_ellipse();

    std::vector<Inference> get_inferences() override;

    void add_inferences_circle_or_ellipse(std::vector<Inference>& inferences);
    void add_inferences_line(std::vector<Inference>& inferences);
    void add_inferences_other(std::vector<Inference>& inferences);

    int _id;
    int _curve;

    void sample_points(
        const int num_points,
        const bool sample_tangents,
        std::vector<Eigen::VectorXd>& samples,
        Eigen::Vector2d& t_range) override;
    
    bool sample_curve(
        const int N_samples, // N
        Eigen::Vector2d& t_bounds,
        Eigen::MatrixXd& t_samples // (Nx6) x,y,z,t_x,t_y,t_z
    );
};

struct OCCTEdge: public Edge {
    OCCTEdge(const TopoDS_Shape& shape, const TopTools_ListOfShape& faces);

    void init_bb();
    void init_nabb();

    void init_line();
    void init_circle();
    void init_ellipse();

    std::vector<Inference> get_inferences() override;

    void add_inferences_circle_or_ellipse(std::vector<Inference>& inferences);
    void add_inferences_line(std::vector<Inference>& inferences);
    void add_inferences_other(std::vector<Inference>& inferences);

    TopoDS_Edge _shape;
    BRepAdaptor_Curve _curve;
    TopTools_ListOfShape _faces;

    void sample_points(
        const int num_points,
        const bool sample_tangents,
        std::vector<Eigen::VectorXd>& samples,
        Eigen::Vector2d& t_range) override;

    bool sample_curve(
        const int N_samples, // N
        Eigen::Vector2d& t_bounds,
        Eigen::MatrixXd& t_samples // (Nx6) x,y,z,t_x,t_y,t_z
    );
};

}

#endif // !EDGE_H_INCLUDED
