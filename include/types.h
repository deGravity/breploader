#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED 1

#include <Eigen/Core>
#include <assert.h>
#include <parasolid.h>
#include <gp_Pnt.hxx>
#include <gp_Mat.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Shape.hxx>
#include <TopAbs_ShapeEnum.hxx>

namespace pspy {

enum class TopologyType {
    FACE,
    EDGE,
    VERTEX,
    LOOP
};

enum class SurfaceFunction {
    PLANE,
    CYLINDER,
    CONE,
    SPHERE,
    TORUS,
    SPUN,
    BSURF,
    OFFSET,
    SWEPT,
    BLENDSF,
    MESH,
    FSURF,
    SURFACEOFEXTRUSION,
    OTHERSURFACE,
    NONE // Faces sometimes have no surface
};

// TODO: clean up potential duplicates introduced by adding OCCT curve functions
enum class CurveFunction {
    LINE,
    CIRCLE,
    ELLIPSE,
    BCURVE,
    ICURVE,
    FCURVE,
    SPCURVE,
    TRCURVE,
    CPCURVE,
    PLINE,
    HYPERBOLA,
    PARABOLA,
    OFFSETCURVE,
    OTHERCURVE,
    NONE // Edges sometimes have no curve
};

enum class LoopType {
    VERTEX,
    WIRE,
    OUTER,
    INNER,
    WINDING,
    INNER_SING,
    LIKELY_OUTER,
    LIKELY_INNER,
    UNCLEAR,
    ERROR
};

enum class InferenceType {
    CENTER, // CIRCLE, ELLIPSE, SPHERE
    CENTROID, // PLANE - ignores inner loops (COG centroid)
    MID_POINT, // EDGES that aren't circles or ellipses
    POINT, // VERTEX, tip of CONE
    TOP_AXIS_POINT, // CYLINDER, CONE, TORUS, SPUN
    BOTTOM_AXIS_POINT, // ...
    MID_AXIS_POINT, // ...
    LOOP_CENTER // First edge of planar loop, also not CIRCLE
};

const double XFRM_TOL = 0.999;
const double MASS_ACC = 0.999;

struct MassProperties {
    MassProperties(int* ids, double accuracy = MASS_ACC, int num_ids = 1) {
        PK_ERROR_t err = PK_ERROR_no_errors;
        PK_TOPOL_eval_mass_props_o_t mass_props_options;
        PK_TOPOL_eval_mass_props_o_m(mass_props_options);

        m_of_i.resize(3, 3);

        err = PK_TOPOL_eval_mass_props(num_ids, ids, accuracy,
            &mass_props_options,
            &amount,
            &mass,
            c_of_g.data(),
            m_of_i.data(),
            &periphery);
        assert(err == PK_ERROR_no_errors || err == PK_ERROR_missing_geom); // PK_TOPOL_eval_mass_props
        if (err == PK_ERROR_missing_geom) {
            amount = 0;
            mass = 0;
            c_of_g.setZero();
            m_of_i = Eigen::MatrixXd::Zero(3, 3);
            periphery = 0;
        }
    }

    MassProperties(const TopoDS_Shape& shape, double accuracy = MASS_ACC) {
        m_of_i.resize(3, 3);

        if (!shape.IsNull()) {
            gp_Pnt c_of_g_pnt;
            gp_Mat m_of_i_mat;
            GProp_GProps linear_props;
            GProp_GProps surface_props;
            GProp_GProps volume_props;

            // TODO: handle different densities.
            //  BRepGProp functions do not attach density to the entities considered,
            //  which may lead to different results from PK_TOPOL_eval_mass_props,
            //  as per
            //  https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html.
            switch (shape.ShapeType()) {
            case TopAbs_COMPOUND:
            case TopAbs_COMPSOLID:
            case TopAbs_SOLID:
            case TopAbs_SHELL:
                BRepGProp::VolumeProperties(shape, volume_props, accuracy);
                BRepGProp::SurfaceProperties(shape, surface_props, accuracy);
                amount = volume_props.Mass();
                periphery = surface_props.Mass();
                c_of_g_pnt = volume_props.CentreOfMass();
                m_of_i_mat = volume_props.MatrixOfInertia();
                break;
            case TopAbs_FACE:
                BRepGProp::SurfaceProperties(shape, surface_props, accuracy);
                BRepGProp::LinearProperties(shape, linear_props);
                amount = surface_props.Mass();
                periphery = linear_props.Mass();
                c_of_g_pnt = surface_props.CentreOfMass();
                m_of_i_mat = surface_props.MatrixOfInertia();
                break;
            case TopAbs_WIRE:
            case TopAbs_EDGE:
                BRepGProp::LinearProperties(shape, linear_props);
                amount = linear_props.Mass();
                periphery = 0;
                c_of_g_pnt = linear_props.CentreOfMass();
                m_of_i_mat = linear_props.MatrixOfInertia();
                break;
            default:
                amount = 0;
                periphery = 0;
                break;
            }

            mass = amount;
            c_of_g << c_of_g_pnt.X(), c_of_g_pnt.Y(), c_of_g_pnt.Z();
            m_of_i <<
                m_of_i_mat.Value(1, 1), m_of_i_mat.Value(1, 2), m_of_i_mat.Value(1, 3),
                m_of_i_mat.Value(2, 1), m_of_i_mat.Value(2, 2), m_of_i_mat.Value(2, 3),
                m_of_i_mat.Value(3, 1), m_of_i_mat.Value(3, 2), m_of_i_mat.Value(3, 3);
        } else {
            amount = 0;
            mass = 0;
            c_of_g.setZero();
            m_of_i = Eigen::MatrixXd::Zero(3, 3);
            periphery = 0;
        }
    }

    double amount;
    double mass;
    Eigen::Vector3d c_of_g;
    Eigen::MatrixXd m_of_i;
    double periphery;
};

struct Inference {
    Eigen::Vector3d z_axis;
    Eigen::Vector3d origin;
    InferenceType inference_type;

    // For interfacing with the Onshape CAD System
    // Onshape has special case rules for some coordinate
    // inferences that can exclude an inference for redundancy
    // (e.g. all circular inner loops), or flip the axis
    // (e.g. circular or elliptical curves oriented opposite
    // of a planar face containing them)
    bool onshape_inference = true; // Onshape allows this inference
    bool flipped_in_onshape = false; // Onshape would flip our inference
};

template <typename T>
struct TopoDS_Shape_Hash {
    size_t operator()(const T& shape) const {
        return shape.HashCode(INT_MAX);
    }
};

template <typename T>
struct TopoDS_Shape_Pred {
    size_t operator()(const T& shape1, const T& shape2) const {
        return shape1.IsSame(shape2);
    }
};

// From https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
struct gp_Pnt_Hash {
    size_t operator()(const gp_Pnt& pnt) const {
        size_t seed = 0;
        for (size_t i = 1; i <= 3; ++i) {
            double elem = pnt.Coord(i);
            seed ^= std::hash<double>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct gp_Pnt_Pred {
    size_t operator()(const gp_Pnt& pnt1, const gp_Pnt& pnt2) const {
        return
            pnt1.X() == pnt2.X() &&
            pnt1.Y() == pnt2.Y() &&
            pnt1.Z() == pnt2.Z();
    }
};

}

#endif // !TYPES_H_INCLUDED