#include <edge.h>

#include <assert.h>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <GProp_GProps.hxx>
#include <GeomLProp_CLProps.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <Bnd_OBB.hxx>
#include <BRepBndLib.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <gp_Pln.hxx>
#include <gp_Circ.hxx>
#include <gp_Lin.hxx>
#include <gp_Elips.hxx>

namespace pspy {

    OCCTEdge::OCCTEdge(const TopoDS_Shape& shape, const TopTools_ListOfShape& faces) {
        assert(shape.ShapeType() == TopAbs_EDGE);

        _shape = TopoDS::Edge(shape);
        _faces = faces;

        gp_Pnt first, last, mid;
        double t_first, t_last, t_mid;
        opencascade::handle<Geom_Curve> curve_handle = BRep_Tool::Curve(_shape, t_first, t_last);

        if (_shape.IsNull() || curve_handle.IsNull()) {
            _curve = BRepAdaptor_Curve();
            _has_curve = false;
            _is_reversed = false;
            start = Eigen::Vector3d(0, 0, 0);
            end = Eigen::Vector3d(0, 0, 0);
            t_start = 0;
            t_end = 0;
            is_periodic = false;
            mid_point = Eigen::Vector3d(0, 0, 0);

            length = 0;
            center_of_gravity.setZero();
            moment_of_inertia = Eigen::MatrixXd::Zero(3, 3);
        }
        else {
            _curve = BRepAdaptor_Curve(_shape);
            _has_curve = true;

            t_mid = (t_first + t_last) / 2;
            curve_handle->D0(t_first, first);
            curve_handle->D0(t_last, last);
            curve_handle->D0(t_mid, mid);

            assert(_shape.Orientation() == TopAbs_FORWARD || _shape.Orientation() == TopAbs_REVERSED);
            if (_shape.Orientation() == TopAbs_FORWARD) {
                start = Eigen::Vector3d(first.X(), first.Y(), first.Z());
                end = Eigen::Vector3d(last.X(), last.Y(), last.Z());
                t_start = t_first;
                t_end = t_last;
            }
            else {
                end = Eigen::Vector3d(first.X(), first.Y(), first.Z());
                start = Eigen::Vector3d(last.X(), last.Y(), last.Z());
                t_end = t_first;
                t_start = t_last;
                _is_reversed = true; // Needed for tangent evaluation
            }

            is_periodic = _curve.IsPeriodic();

            mid_point = Eigen::Vector3d(mid.X(), mid.Y(), mid.Z());

            auto m = MassProperties(_shape);
            length = m.amount;
            center_of_gravity = m.c_of_g;
            moment_of_inertia = m.m_of_i;
        }

        init_bb();

        init_nabb();

        if (!_has_curve) {
            function = CurveFunction::NONE;
            return;
        }

        GeomAbs_CurveType curve_class = _curve.GetType();

        switch (curve_class) {
        case GeomAbs_Line:
            init_line();
            break;
        case GeomAbs_Circle:
            init_circle();
            break;
        case GeomAbs_Ellipse:
            init_ellipse();
            break;
        case GeomAbs_Hyperbola:
            function = CurveFunction::HYPERBOLA;
            break;
        case GeomAbs_Parabola:
            function = CurveFunction::PARABOLA;
            break;
        case GeomAbs_BezierCurve:
            function = CurveFunction::BCURVE;
            break;
        case GeomAbs_BSplineCurve:
            function = CurveFunction::BCURVE;
            break;
        case GeomAbs_OffsetCurve:
            function = CurveFunction::OFFSETCURVE;
            break;
        case GeomAbs_OtherCurve:
            function = CurveFunction::OTHERCURVE;
            break;
        }

    }

    void OCCTEdge::init_line() {
        gp_Lin line = _curve.Line();
        function = CurveFunction::LINE;
        parameters.push_back(line.Location().X());
        parameters.push_back(line.Location().Y());
        parameters.push_back(line.Location().Z());
        parameters.push_back(line.Direction().X());
        parameters.push_back(line.Direction().Y());
        parameters.push_back(line.Direction().Z());
    }

    void OCCTEdge::init_circle() {
        gp_Circ circle = _curve.Circle();
        function = CurveFunction::CIRCLE;
        parameters.push_back(circle.Location().X());
        parameters.push_back(circle.Location().Y());
        parameters.push_back(circle.Location().Z());
        parameters.push_back(circle.Axis().Direction().X());
        parameters.push_back(circle.Axis().Direction().Y());
        parameters.push_back(circle.Axis().Direction().Z());
        parameters.push_back(circle.XAxis().Direction().X());
        parameters.push_back(circle.XAxis().Direction().Y());
        parameters.push_back(circle.XAxis().Direction().Z());
        parameters.push_back(circle.Radius());
    }

    void OCCTEdge::init_ellipse() {
        gp_Elips ellipse = _curve.Ellipse();
        function = CurveFunction::ELLIPSE;
        parameters.push_back(ellipse.Location().X());
        parameters.push_back(ellipse.Location().Y());
        parameters.push_back(ellipse.Location().Z());
        parameters.push_back(ellipse.Axis().Direction().X());
        parameters.push_back(ellipse.Axis().Direction().Y());
        parameters.push_back(ellipse.Axis().Direction().Z());
        parameters.push_back(ellipse.XAxis().Direction().X());
        parameters.push_back(ellipse.XAxis().Direction().Y());
        parameters.push_back(ellipse.XAxis().Direction().Z());
        parameters.push_back(ellipse.MinorRadius());
        parameters.push_back(ellipse.MajorRadius());
    }

    void OCCTEdge::init_bb() {
        if (!_has_curve) {
            bounding_box = Eigen::MatrixXd::Zero(2, 3);
        }
        else {
            Bnd_Box box;
            BRepBndLib::Add(_shape, box);

            double xmin, ymin, zmin, xmax, ymax, zmax;
            box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

            bounding_box.resize(2, 3);
            bounding_box <<
                xmin, ymin, zmin,
                xmax, ymax, zmax;
        }
    }

    void OCCTEdge::init_nabb() {
        // Get Non-Aligned (Oriented) Bounding Box
        // Axes are ordered so X is largest, Y is second largest, Z is smallest

        if (!_has_curve) {
            na_bb_center.setZero();
            na_bb_x.setZero();
            na_bb_z.setZero();
            na_bounding_box = Eigen::MatrixXd::Zero(2, 3);
            return;
        }
        else {
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
    }

    std::vector<Inference> OCCTEdge::get_inferences()
    {
        std::vector<Inference> inferences;

        switch (function) {
        case CurveFunction::CIRCLE:
            add_inferences_circle_or_ellipse(inferences);
            break;
        case CurveFunction::ELLIPSE:
            add_inferences_circle_or_ellipse(inferences);
            break;
        case CurveFunction::LINE:
            add_inferences_line(inferences);
            break;
        default:
            add_inferences_other(inferences);
            break;
        }

        return inferences;
    }

    void OCCTEdge::add_inferences_circle_or_ellipse(std::vector<Inference>& inferences)
    {
        Inference inf;
        inf.inference_type = InferenceType::CENTER;
        inf.origin = Eigen::Vector3d(parameters[0], parameters[1], parameters[2]);
        inf.z_axis = Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);

        // Onshape Special Case: If a circle or ellipse is on a planar face,
        // orient with the face normal instead of the cirlce or ellipse axis
        auto is_parallel = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
            double dot_product = a.dot(b);
            double diff = (1 - (dot_product * dot_product));
            double abs_diff = diff * diff;
            return (abs_diff < FLT_EPSILON);
        };

        for (auto iterator = _faces.cbegin(); iterator != _faces.cend(); iterator++) {
            TopoDS_Face face_shape = TopoDS::Face(*iterator);
            BRepAdaptor_Surface surf(face_shape);
            TopAbs_Orientation orientation = face_shape.Orientation();
            GeomAbs_SurfaceType face_class = surf.GetType();

            if (face_class == GeomAbs_Plane) {
                gp_Pln plane = surf.Plane();
                Eigen::Vector3d plane_axis(
                    plane.Axis().Direction().X(),
                    plane.Axis().Direction().Y(),
                    plane.Axis().Direction().Z());
                if (orientation == TopAbs_REVERSED) {
                    plane_axis = -plane_axis;
                }

                // If parallel then the circle lies in the face
                if (is_parallel(plane_axis, inf.z_axis)) {
                    // If parallel and not equal, then they are opposite
                    // and we should flip the axis if using with Onshape
                    if ((plane_axis.normalized() -
                        inf.z_axis.normalized()).norm() > FLT_EPSILON) {
                        inf.flipped_in_onshape = true;
                    }
                    break; // Stop after we find the containing plane, if it exists
                }
            }
        }

        inferences.push_back(inf);
    }

    void OCCTEdge::add_inferences_line(std::vector<Inference>& inferences)
    {
        Inference inf;
        inf.inference_type = InferenceType::MID_POINT;
        inf.origin = mid_point;
        inf.z_axis = Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);

        inferences.push_back(inf);
    }

    void OCCTEdge::add_inferences_other(std::vector<Inference>& inferences)
    {
        Inference inf;
        inf.inference_type = InferenceType::MID_POINT;
        inf.origin = mid_point;
        inf.z_axis = Eigen::Vector3d(0.0, 0.0, 1.0);

        inferences.push_back(inf);
    }

    void OCCTEdge::sample_points(
        const int num_points,
        const bool sample_tangents,
        std::vector<Eigen::VectorXd>& samples,
        Eigen::Vector2d& t_range) {
        // TODO: implement
    }

    bool OCCTEdge::sample_curve(
        const int N_samples,
        Eigen::Vector2d& t_bounds,
        Eigen::MatrixXd& t_samples
    ) {
        // Start and end points on the _edge_ may be reversed from the original curve
        t_bounds(0) = _is_reversed ? t_end : t_start;
        t_bounds(1) = _is_reversed ? t_start : t_end;

        // Only try evaluating a curve if it exists
        if (!_has_curve) {
            return false;
        }

        double t_first, t_last;
        opencascade::handle<Geom_Curve> curve_handle = BRep_Tool::Curve(_shape, t_first, t_last);

        Eigen::ArrayXd ts = Eigen::ArrayXd::LinSpaced(N_samples, t_bounds(0), t_bounds(1));
        t_samples.resize(N_samples, 6);
        gp_Pnt point;
        gp_Vec tangent;
        gp_Dir tangent_dir;
        for (int i = 0; i < N_samples; ++i) {
            // Sample the curve at ts(i), with no derivatives
            GeomLProp_CLProps curve_props(curve_handle, ts(i), 0, DBL_EPSILON);
            point = curve_props.Value();

            if (!curve_props.IsTangentDefined()) {
                // Use (0,0,0) as singularity tangent
                tangent.SetCoord(1, 0);
                tangent.SetCoord(2, 0);
                tangent.SetCoord(3, 0);
            }
            else {
                // Use the actual tangent
                curve_props.Tangent(tangent_dir);
                tangent = gp_Vec(tangent_dir);
            }
            for (int j = 0; j < 3; ++j) {
                t_samples(i, j) = point.Coord(j + 1);
                t_samples(i, j + 3) = tangent.Coord(j + 1);
            }
        }
        return true;
    }

}