#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED 1

#include <vector>
#include <Eigen/Core>
#include <TopoDS_Face.hxx>
#include <BRepAdaptor_Surface.hxx>
#include "types.h"

namespace pspy {

struct Face {
    virtual std::vector<Inference> get_inferences() = 0;

    bool _has_surface;
    SurfaceFunction function;
    std::vector<double> parameters;
    bool orientation;

    Eigen::MatrixXd bounding_box;
    Eigen::Vector3d na_bb_center;
    Eigen::Vector3d na_bb_x;
    Eigen::Vector3d na_bb_z;
    Eigen::MatrixXd na_bounding_box;

    double surface_area;
    double circumference;

    Eigen::Vector3d center_of_gravity;
    Eigen::MatrixXd moment_of_inertia;

    virtual void sample_points(
        const int num_points,
        const bool sample_normals,
        std::vector<Eigen::MatrixXd>& samples,
        Eigen::MatrixXd& uv_box) = 0;

    virtual void random_sample_points(
        const int num_points,
        Eigen::MatrixXd& samples,
        Eigen::MatrixXd& coords,
        Eigen::MatrixXd& uv_box) = 0;
    
    virtual bool sample_surface(
        const int N_ref_samples,
        const int N_uv_samples, // N
        Eigen::MatrixXd& uv_bounds,
        Eigen::MatrixXd& uv_coords, // (Nx2) Sampled Coords normalized to [0,1]
        Eigen::MatrixXd& uv_samples // (Nx7) x,y,z,n_x,n_y,n_z,m
    ) = 0;

};

struct PSFace: public Face {
    PSFace(int id);

    void init_parametric_function();

    void init_bb();
    void init_nabb();
    void init_mass_props();

    void init_plane();
    void init_cyl();
    void init_cone();
    void init_sphere();
    void init_torus();
    void init_spun();

    std::vector<Inference> get_inferences() override;

    void add_inferences_plane(std::vector<Inference>& inferences);
    void add_inferences_cone(std::vector<Inference>& inferences);
    void add_inferences_sphere(std::vector<Inference>& inferences);
    void add_inferences_axial(std::vector<Inference>& inferences);

    int _id;
    int _surf;

    void sample_points(
        const int num_points,
        const bool sample_normals,
        std::vector<Eigen::MatrixXd>& samples,
        Eigen::MatrixXd& uv_box) override;

    void random_sample_points(
        const int num_points,
        Eigen::MatrixXd& samples,
        Eigen::MatrixXd& coords,
        Eigen::MatrixXd& uv_box) override;

    void sample_mask_sdf(
        const int quality,
        const int num_points,
        Eigen::MatrixXd& coords,
        Eigen::VectorXd& sdf,
        Eigen::MatrixXd& uv_box
    );

    bool sample_surface(
        const int N_ref_samples,
        const int N_uv_samples, // N
        Eigen::MatrixXd& uv_bounds,
        Eigen::MatrixXd& uv_coords, // (Nx2) Sampled Coords normalized to [0,1]
        Eigen::MatrixXd& uv_samples // (Nx7) x,y,z,n_x,n_y,n_z,m
    ) override;
};

struct OCCTFace: public Face {
    OCCTFace(const TopoDS_Shape& shape);

    void init_parametric_function();

    void init_bb();
    void init_nabb();
    void init_mass_props();

    void init_plane();
    void init_cyl();
    void init_cone();
    void init_sphere();
    void init_torus();
    void init_spun();

    std::vector<Inference> get_inferences() override;

    void add_inferences_plane(std::vector<Inference>& inferences);
    void add_inferences_cone(std::vector<Inference>& inferences);
    void add_inferences_sphere(std::vector<Inference>& inferences);
    void add_inferences_axial(std::vector<Inference>& inferences);

    TopoDS_Face _shape;
    BRepAdaptor_Surface _surf;

    void sample_points(
        const int num_points,
        const bool sample_normals,
        std::vector<Eigen::MatrixXd>& samples,
        Eigen::MatrixXd& uv_box) override;

    void random_sample_points(
        const int num_points,
        Eigen::MatrixXd& samples,
        Eigen::MatrixXd& coords,
        Eigen::MatrixXd& uv_box) override;
    
    bool sample_surface(
        const int N_ref_samples,
        const int N_uv_samples, // N
        Eigen::MatrixXd& uv_bounds,
        Eigen::MatrixXd& uv_coords, // (Nx2) Sampled Coords normalized to [0,1]
        Eigen::MatrixXd& uv_samples // (Nx7) x,y,z,n_x,n_y,n_z,m
    ) override;

};

}

#endif // !FACE_H_INCLUDED
