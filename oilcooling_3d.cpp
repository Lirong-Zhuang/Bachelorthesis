/**
 * @file 	oil_cooling_3D.cpp
 * @brief 	3D model oil cooling system in a running electric motor.
 * @author 	Lirong Zhuang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
std::string hausing_CAD = "./input/Electricmotor_Hausing.STL";
std::string rotor_CAD = "./input/Electricmotor_Rotor.STL";
std::string winding_CAD = "./input/Electricmotor_Winding.STL";
std::string cover_CAD = "./input/Electricmotor_Cover.STL";
// std::string fluid_CAD = "./input/Electricmotor_Fluid.STL";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.25;               // Domain length.
Real DH = 0.25;               // Domain height.
Real DW = 0.25;               // Domain width.
Real DO = 0.02;               /**< Outflow diameter. */
Real resolution_ref = 0.0028; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4; /**< Extending width for wall boundary. */
Vec3d domain_lower_bound(-DL, -DH, -DW);
Vec3d domain_upper_bound(DL, DH, DW);
//Vecd translation(-0.129, -0.127, -0.12);
Vecd translation(-0.11, -0.11, -0.12);
Real scaling = 0.001;
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Geometrie parameters.
//----------------------------------------------------------------------
Real Wnum = 12; /**< Winding Number. */
Real angle_increment = 2 * Pi / Wnum;
Real DMO = 0.22;                     /**< Diameter of outer motor hausing. */
Real RMO = 0.5 * DMO;                /**< Radius of outer motor hausing. */
Real DMN = 0.196;                    /**< Diameter of inner motor hausing. */
Real RMI = 0.5 * DMN;                /**< Radius of inner motor hausing. */
Real DRO = 0.12;                     /**< Diameter of outer rotor. */
Real LL = 0.006;                     /**< Inflow region length. */
Real LH = RMO - RMI;                 /**< Inflows region height. */
Real WL = 0.03;                      /**< Winding length. */
Real WH = 0.025;                     /**< Winding height. */
Real WLK = 0.009;                    /**< Winding thickness. */
Real WCL = WL - 2 * WLK;             /**< Winding core thickness. */
Real AG = 0.001;                     /**< Air-gap. */
Real WD = 0.5 * DRO + AG + 0.5 * WH; /**< Distance from center point to the center of a winding. */
Real inlet_height = DMO - LH;        /**< Inflow location height */
Real rho0_f = 812;                   /**< Reference density of fluid. */
Real gravity_g = 9.81;               /**< Gravity force of fluid. */
// Axis
Vec3d axis_x(1.0, 0.0, 0.0);
Vec3d axis_y(0.0, 1.0, 0.0);
Vec3d axis_z(0.0, 0.0, 1.0);
// dynamics informations of oil
Real flow_rate = 110;                                                     /**< Oil flow rate [L / h] */
Real v_inlet = (flow_rate * 0.004 / 3600) / (Pi * LL * LL);               /**< Inflow vilocity [m / s]. */
Real U_f = sqrt(v_inlet * v_inlet + 2 * gravity_g * (inlet_height + LH)); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                                    /**< Reference sound speed. */
// dynamics informations of rotor
Real rotor_rotation_velocity = 900;                    /**<Angular velocity rpm. */
Real Omega = -(rotor_rotation_velocity * 2 * Pi / 60); /**<Angle of rotor. */
// thermal parameters
Real mu_f = 0.00975;                                            /**< Dynamics viscosity [Pa * s]. */
Real phi_winding = 110.0;                                       /**< Temperature of winding at begin. */
Real phi_fluid_initial = 75.0;                                  /**< Temperature of oil at begin. */
Real k_oil = 7.62e-8;                                           /**< Diffusion coefficient of oil 2.0e-7. */
Real k_winding = 1.14e-4;                                       /**< Diffusion coefficient of winding 1.1e-6. */
Real k_contact = (2 * k_oil * k_winding) / (k_oil + k_winding); /**< Thermal conductivity between winding and oil. */
Real dq = 0.5;                                                  /**< Heating efficient of internal heat source [°C/s]. */
//----------------------------------------------------------------------
//	Geometrie of the inlets and oulets.
//----------------------------------------------------------------------
Real R_inlet = RMO - 0.5 * LH;
Real L_inlet = 0.0745;
Real R_outlet = -(0.5 * (RMO + RMI));
Real L_outlet = 0.085;
Vec3d inlet_halfsize(0.5 * LH, 0.5 * LL, 0.5 * LL);
Vec3d outlet_halfsize(0.5 * (RMO - RMI), 0.5 * DO, 0.5 * DO);
Vec3d inleta1_translation = Vec3d(0, R_inlet, L_inlet);
Transform inleta1_transform(Rotation3d(-(Pi / 2), axis_z), inleta1_translation);
Vec3d inleta2_translation = Vec3d(R_inlet * cos(5 * angle_increment), R_inlet *sin(5 * angle_increment), L_inlet);
Transform inleta2_transform(Rotation3d(-(angle_increment), axis_z), inleta2_translation);
Vec3d inleta3_translation = Vec3d(R_inlet * cos(4 * angle_increment), R_inlet *sin(4 * angle_increment), L_inlet);
Transform inleta3_transform(Rotation3d(-(2 * angle_increment), axis_z), inleta3_translation);
Vec3d inleta4_translation = Vec3d(R_inlet * cos(2 * angle_increment), R_inlet *sin(2 * angle_increment), L_inlet);
Transform inleta4_transform(Rotation3d(-(4 * angle_increment), axis_z), inleta4_translation);
Vec3d inleta5_translation = Vec3d(R_inlet * cos(angle_increment), R_inlet *sin(angle_increment), L_inlet);
Transform inleta5_transform(Rotation3d(-(5 * angle_increment), axis_z), inleta5_translation);
Vec3d outleta_translation = Vec3d(0, R_outlet, L_outlet);
Transform outleta_transform(Rotation3d(-(Pi / 2), axis_z), outleta_translation);
Vec3d inletb1_translation = Vec3d(0, R_inlet, -(L_inlet));
Transform inletb1_transform(Rotation3d(-(Pi / 2), axis_z), inletb1_translation);
Vec3d inletb2_translation = Vec3d(R_inlet * cos(5 * angle_increment), R_inlet *sin(5 * angle_increment), -(L_inlet));
Transform inletb2_transform(Rotation3d(-(angle_increment), axis_z), inletb2_translation);
Vec3d inletb3_translation = Vec3d(R_inlet * cos(4 * angle_increment), R_inlet *sin(4 * angle_increment), -(L_inlet));
Transform inletb3_transform(Rotation3d(-(2 * angle_increment), axis_z), inletb3_translation);
Vec3d inletb4_translation = Vec3d(R_inlet * cos(2 * angle_increment), R_inlet *sin(2 * angle_increment), -(L_inlet));
Transform inletb4_transform(Rotation3d(-(4 * angle_increment), axis_z), inletb4_translation);
Vec3d inletb5_translation = Vec3d(R_inlet * cos(angle_increment), R_inlet *sin(angle_increment), -(L_inlet));
Transform inletb5_transform(Rotation3d(-(5 * angle_increment), axis_z), inletb5_translation);
Vec3d outletb_translation = Vec3d(0, R_outlet, -(L_outlet));
Transform outletb_transform(Rotation3d(-(Pi / 2), axis_z), outletb_translation);
SimTK::UnitVec3 axis1(0.0, 1.0, 0.0);
SimTK::UnitVec3 axis2(cos(5 * angle_increment), sin(5 * angle_increment), 0.0);
SimTK::UnitVec3 axis3(cos(4 * angle_increment), sin(4 * angle_increment), 0.0);
SimTK::UnitVec3 axis4(cos(2 * angle_increment), sin(2 * angle_increment), 0.0);
SimTK::UnitVec3 axis5(cos(1 * angle_increment), sin(1 * angle_increment), 0.0);
//----------------------------------------------------------------------
//	Locations of observer.
//----------------------------------------------------------------------
Real value_z_slots = 0.07;   /**< Location of slots at axis z(if all 12 slots are in the same xy-flat. */
Real value_z_toothes = 0.07; /**< Location of slots at axis z(if all 12 toothes are in the same xy-flat. */
Real value_z_yokes = 0.07;   /**< Location of slots at axis z(if all 12 yokes are in the same xy-flat. */
StdVec<Vec3d> slot_locations = {
    Vec3d(WD, 0.5 * WCL, value_z_slots),
    Vec3d(WD *cos(angle_increment) - 0.5 * WCL * sin(angle_increment),
          WD *sin(angle_increment) + 0.5 * WCL * cos(angle_increment),
          value_z_slots),
    Vec3d(WD *cos(2 * angle_increment) - 0.5 * WCL * sin(2 * angle_increment),
          WD *sin(2 * angle_increment) + 0.5 * WCL * cos(2 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(3 * angle_increment) - 0.5 * WCL * sin(3 * angle_increment),
          WD *sin(3 * angle_increment) + 0.5 * WCL * cos(3 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(4 * angle_increment) - 0.5 * WCL * sin(4 * angle_increment),
          WD *sin(4 * angle_increment) + 0.5 * WCL * cos(4 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(5 * angle_increment) - 0.5 * WCL * sin(5 * angle_increment),
          WD *sin(5 * angle_increment) + 0.5 * WCL * cos(5 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(6 * angle_increment) - 0.5 * WCL * sin(6 * angle_increment),
          WD *sin(6 * angle_increment) + 0.5 * WCL * cos(6 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(7 * angle_increment) - 0.5 * WCL * sin(7 * angle_increment),
          WD *sin(7 * angle_increment) + 0.5 * WCL * cos(7 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(8 * angle_increment) - 0.5 * WCL * sin(8 * angle_increment),
          WD *sin(8 * angle_increment) + 0.5 * WCL * cos(8 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(9 * angle_increment) - 0.5 * WCL * sin(9 * angle_increment),
          WD *sin(9 * angle_increment) + 0.5 * WCL * cos(9 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(10 * angle_increment) - 0.5 * WCL * sin(10 * angle_increment),
          WD *sin(10 * angle_increment) + 0.5 * WCL * cos(10 * angle_increment),
          value_z_slots),
    Vec3d(WD *cos(11 * angle_increment) - 0.5 * WCL * sin(11 * angle_increment),
          WD *sin(11 * angle_increment) + 0.5 * WCL * cos(11 * angle_increment),
          value_z_slots)
};
StdVec<Vec3d> tooth_locations = {
    Vec3d(WD, 0.5 * WL, value_z_toothes),
    Vec3d(WD *cos(angle_increment) - 0.5 * WL * sin(angle_increment),
          WD *sin(angle_increment) + 0.5 * WL * cos(angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(2 * angle_increment) - 0.5 * WL * sin(2 * angle_increment),
          WD *sin(2 * angle_increment) + 0.5 * WL * cos(2 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(3 * angle_increment) - 0.5 * WL * sin(3 * angle_increment),
          WD *sin(3 * angle_increment) + 0.5 * WL * cos(3 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(4 * angle_increment) - 0.5 * WL * sin(4 * angle_increment),
          WD *sin(4 * angle_increment) + 0.5 * WL * cos(4 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(5 * angle_increment) - 0.5 * WL * sin(5 * angle_increment),
          WD *sin(5 * angle_increment) + 0.5 * WL * cos(5 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(6 * angle_increment) - 0.5 * WL * sin(6 * angle_increment),
          WD *sin(6 * angle_increment) + 0.5 * WL * cos(6 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(7 * angle_increment) - 0.5 * WL * sin(7 * angle_increment),
          WD *sin(7 * angle_increment) + 0.5 * WL * cos(7 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(8 * angle_increment) - 0.5 * WL * sin(8 * angle_increment),
          WD *sin(8 * angle_increment) + 0.5 * WL * cos(8 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(9 * angle_increment) - 0.5 * WL * sin(9 * angle_increment),
          WD *sin(9 * angle_increment) + 0.5 * WL * cos(9 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(10 * angle_increment) - 0.5 * WL * sin(10 * angle_increment),
          WD *sin(10 * angle_increment) + 0.5 * WL * cos(10 * angle_increment),
          value_z_toothes),
    Vec3d(WD *cos(11 * angle_increment) - 0.5 * WL * sin(11 * angle_increment),
          WD *sin(11 * angle_increment) + 0.5 * WL * cos(11 * angle_increment),
          value_z_toothes)
};
StdVec<Vec3d> yoke_locations = {
    Vec3d(WD + 0.5 * WH, 0.5 * WCL + 0.5 * WLK, value_z_yokes),
    Vec3d(WD *cos(angle_increment) + 0.5 * WH * cos(angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(angle_increment),
          WD *sin(angle_increment) + 0.5 * WH * sin(angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(2 * angle_increment) + 0.5 * WH * cos(2 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(2 * angle_increment),
          WD *sin(2 * angle_increment) + 0.5 * WH * sin(2 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(2 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(3 * angle_increment) + 0.5 * WH * cos(3 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(3 * angle_increment),
          WD *sin(3 * angle_increment) + 0.5 * WH * sin(3 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(3 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(4 * angle_increment) + 0.5 * WH * cos(4 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(4 * angle_increment),
          WD *sin(4 * angle_increment) + 0.5 * WH * sin(4 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(4 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(5 * angle_increment) + 0.5 * WH * cos(5 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(5 * angle_increment),
          WD *sin(5 * angle_increment) + 0.5 * WH * sin(5 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(5 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(6 * angle_increment) + 0.5 * WH * cos(6 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(6 * angle_increment),
          WD *sin(6 * angle_increment) + 0.5 * WH * sin(6 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(6 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(7 * angle_increment) + 0.5 * WH * cos(7 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(7 * angle_increment),
          WD *sin(7 * angle_increment) + 0.5 * WH * sin(7 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(7 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(8 * angle_increment) + 0.5 * WH * cos(8 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(8 * angle_increment),
          WD *sin(8 * angle_increment) + 0.5 * WH * sin(8 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(8 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(9 * angle_increment) + 0.5 * WH * cos(9 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(9 * angle_increment),
          WD *sin(9 * angle_increment) + 0.5 * WH * sin(9 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(9 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(10 * angle_increment) + 0.5 * WH * cos(10 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(10 * angle_increment),
          WD *sin(10 * angle_increment) + 0.5 * WH * sin(10 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(10 * angle_increment),
          value_z_yokes),
    Vec3d(WD *cos(11 * angle_increment) + 0.5 * WH * cos(11 * angle_increment) - (0.5 * WCL + 0.5 * WLK) * sin(11 * angle_increment),
          WD *sin(11 * angle_increment) + 0.5 * WH * sin(11 * angle_increment) + (0.5 * WCL + 0.5 * WLK) * cos(11 * angle_increment),
          value_z_yokes)
};
//----------------------------------------------------------------------
//	Case-dependent Hausing
//----------------------------------------------------------------------
class HausingFromMesh : public ComplexShape
{
  public:
    explicit HausingFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(0.5 * resolution_ref, hausing_CAD, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Case-dependent Rotor
// -------------------------------------------
class RotorFromMesh : public ComplexShape
{
  public:
    explicit RotorFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(2 * resolution_ref, rotor_CAD, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Case-dependent Winding
//----------------------------------------------------------------------
class WindingFromMesh : public ComplexShape
{
  public:
    explicit WindingFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(resolution_ref, winding_CAD, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Case-dependent Winding
//----------------------------------------------------------------------
class CoverFromMesh : public ComplexShape
{
  public:
    explicit CoverFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(2 * resolution_ref, cover_CAD, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Case-dependent Fluid
//----------------------------------------------------------------------
class FluidFromMesh : public ComplexShape
{
  public:
    explicit FluidFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(axis1, 0.5 * LL, 0.5 * LH, resolution_ref, inleta1_translation);
        add<TriangleMeshShapeCylinder>(axis2, 0.5 * LL, 0.5 * LH, resolution_ref, inleta2_translation);
        add<TriangleMeshShapeCylinder>(axis3, 0.5 * LL, 0.5 * LH, resolution_ref, inleta3_translation);
        add<TriangleMeshShapeCylinder>(axis4, 0.5 * LL, 0.5 * LH, resolution_ref, inleta4_translation);
        add<TriangleMeshShapeCylinder>(axis5, 0.5 * LL, 0.5 * LH, resolution_ref, inleta5_translation);
        add<TriangleMeshShapeCylinder>(axis1, 0.5 * LL, 0.5 * LH, resolution_ref, inletb1_translation);
        add<TriangleMeshShapeCylinder>(axis2, 0.5 * LL, 0.5 * LH, resolution_ref, inletb2_translation);
        add<TriangleMeshShapeCylinder>(axis3, 0.5 * LL, 0.5 * LH, resolution_ref, inletb3_translation);
        add<TriangleMeshShapeCylinder>(axis4, 0.5 * LL, 0.5 * LH, resolution_ref, inletb4_translation);
        add<TriangleMeshShapeCylinder>(axis5, 0.5 * LL, 0.5 * LH, resolution_ref, inletb5_translation);
        // add<ExtrudeShape<TriangleMeshShapeSTL>>(resolution_ref, fluid_CAD, translation, scaling);
    }
};
//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
  public:
    InletInflowCondition(BodyAlignedBoxByParticle &aligned_box_part)
        : EmitterInflowCondition(aligned_box_part) {}

  protected:
    virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
    {
        return Vec3d(v_inlet, 0.0, 0.0);
    }
};
//----------------------------------------------------------------------
//	Application dependent initial condition of winding.
//----------------------------------------------------------------------
class ThermoWindingInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermoWindingInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};
    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = phi_winding;
    };

  protected:
    StdLargeVec<Real> &phi_;
};
//----------------------------------------------------------------------
//	The windings heat up due to the internal heat sources.
//----------------------------------------------------------------------
class ThermoWindingHeatSource : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermoWindingHeatSource(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};
    void update(size_t index_i, Real dt)
    {
        phi_[index_i] += dq * dt;
    }

  protected:
    StdLargeVec<Real> &phi_;
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermofluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = phi_fluid_initial;
    };

  protected:
    StdLargeVec<Real> &phi_;
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
using ThermalRelaxationComplex = DiffusionBodyRelaxationComplex<
    IsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Build up a SPHSystem */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody hausing(sph_system, makeShared<HausingFromMesh>("MotorHausing"));
    hausing.defineMaterial<Solid>();
    hausing.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? hausing.generateParticles<BaseParticles, Reload>(hausing.getName())
        : hausing.generateParticles<BaseParticles, Lattice>();

    SolidBody rotor(sph_system, makeShared<RotorFromMesh>("MotorRotor"));
    rotor.defineMaterial<Solid>();
    rotor.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? rotor.generateParticles<BaseParticles, Reload>(rotor.getName())
        : rotor.generateParticles<BaseParticles, Lattice>();

    SolidBody winding(sph_system, makeShared<WindingFromMesh>("MotorWinding"));
    winding.defineMaterial<Solid>();
    winding.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? winding.generateParticles<BaseParticles, Reload>(winding.getName())
        : winding.generateParticles<BaseParticles, Lattice>();

    SolidBody cover(sph_system, makeShared<CoverFromMesh>("MotorCover"));
    cover.defineMaterial<Solid>();
    cover.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cover.generateParticles<BaseParticles, Reload>(cover.getName())
        : cover.generateParticles<BaseParticles, Lattice>();

    FluidBody oil(sph_system, makeShared<FluidFromMesh>("MotorOil"));
    oil.sph_adaptation_->resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    oil.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(3000.0);
    oil.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    ObserverBody slot_observer(sph_system, "ObserverSlot");
    slot_observer.generateParticles<ObserverParticles>(slot_locations);
    ObserverBody tooth_observer(sph_system, "ObserverTooth");
    tooth_observer.generateParticles<ObserverParticles>(tooth_locations);
    ObserverBody yoke_observer(sph_system, "ObserverYoke");
    yoke_observer.generateParticles<ObserverParticles>(yoke_locations);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation hausing_inner(hausing);
        InnerRelation rotor_inner(rotor);
        InnerRelation winding_inner(winding);
        InnerRelation cover_inner(cover);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_hausing_particles(hausing);
        SimpleDynamics<RandomizeParticlePosition> random_rotor_particles(rotor);
        SimpleDynamics<RandomizeParticlePosition> random_winding_particles(winding);
        SimpleDynamics<RandomizeParticlePosition> random_cover_particles(cover);
        RelaxationStepInner relaxation_step_inner_hausing(hausing_inner);
        RelaxationStepInner relaxation_step_inner_rotor(rotor_inner);
        RelaxationStepInner relaxation_step_inner_winding(winding_inner);
        RelaxationStepInner relaxation_step_inner_cover(cover_inner);
        BodyStatesRecordingToVtp write_hausing_to_vtp(hausing);
        BodyStatesRecordingToVtp write_rotor_to_vtp(rotor);
        BodyStatesRecordingToVtp write_winding_to_vtp(winding);
        BodyStatesRecordingToVtp write_cover_to_vtp(cover);
        ReloadParticleIO write_hausing_particle_reload_files(hausing);
        ReloadParticleIO write_rotor_particle_reload_files(rotor);
        ReloadParticleIO write_winding_particle_reload_files(winding);
        ReloadParticleIO write_cover_particle_reload_files(cover);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_hausing_particles.exec(0.25);
        random_rotor_particles.exec(0.25);
        random_winding_particles.exec(0.25);
        random_cover_particles.exec(0.25);
        relaxation_step_inner_hausing.SurfaceBounding().exec();
        relaxation_step_inner_rotor.SurfaceBounding().exec();
        relaxation_step_inner_winding.SurfaceBounding().exec();
        relaxation_step_inner_cover.SurfaceBounding().exec();
        write_hausing_to_vtp.writeToFile(0);
        write_rotor_to_vtp.writeToFile(0);
        write_winding_to_vtp.writeToFile(0);
        write_cover_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner_hausing.exec();
            relaxation_step_inner_rotor.exec();
            relaxation_step_inner_winding.exec();
            relaxation_step_inner_cover.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_hausing_to_vtp.writeToFile(ite_p);
                write_rotor_to_vtp.writeToFile(ite_p);
                write_winding_to_vtp.writeToFile(ite_p);
                write_cover_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_hausing_particle_reload_files.writeToFile(0);
        write_rotor_particle_reload_files.writeToFile(0);
        write_winding_particle_reload_files.writeToFile(0);
        write_cover_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation oil_inner(oil);
    InnerRelation winding_inner(winding);
    ContactRelation oil_solid_contact(oil, {&hausing, &rotor, &winding, &cover});
    ContactRelation oil_thermal_contact(oil, {&winding});
    ContactRelation winding_thermal_contact(winding, {&oil});
    ContactRelation winding_slot_contact(slot_observer, {&winding});
    ContactRelation winding_tooth_contact(tooth_observer, {&winding});
    ContactRelation winding_yoke_contact(yoke_observer, {&winding});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation oil_body_complex(oil_inner, oil_solid_contact);
    ComplexRelation oil_thermal_complex(oil_inner, oil_thermal_contact);
    ComplexRelation winding_thermal_complex(winding_inner, winding_thermal_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> hausing_normal_direction(hausing);
    SimpleDynamics<NormalDirectionFromBodyShape> rotor_normal_direction(rotor);
    SimpleDynamics<NormalDirectionFromBodyShape> winding_normal_direction(winding);
    SimpleDynamics<NormalDirectionFromBodyShape> cover_normal_direction(cover);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(oil_inner, oil_solid_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(oil_inner, oil_solid_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(oil_inner, oil_solid_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(oil_inner, oil_solid_contact);

    Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity(oil, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(oil, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(oil);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(oil_inner, oil_solid_contact);

    BodyAlignedBoxByParticle emitter_a1(oil, makeShared<AlignedBoxShape>(xAxis, inleta1_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_a1(emitter_a1);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_a1(emitter_a1, inlet_buffer);
    BodyAlignedBoxByParticle emitter_a2(oil, makeShared<AlignedBoxShape>(xAxis, inleta2_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_a2(emitter_a2);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_a2(emitter_a2, inlet_buffer);
    BodyAlignedBoxByParticle emitter_a3(oil, makeShared<AlignedBoxShape>(xAxis, inleta3_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_a3(emitter_a3);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_a3(emitter_a3, inlet_buffer);
    BodyAlignedBoxByParticle emitter_a4(oil, makeShared<AlignedBoxShape>(xAxis, inleta4_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_a4(emitter_a4);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_a4(emitter_a4, inlet_buffer);
    BodyAlignedBoxByParticle emitter_a5(oil, makeShared<AlignedBoxShape>(xAxis, inleta5_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_a5(emitter_a5);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_a5(emitter_a5, inlet_buffer);
    BodyAlignedBoxByParticle emitter_b1(oil, makeShared<AlignedBoxShape>(xAxis, inletb1_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_b1(emitter_b1);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_b1(emitter_b1, inlet_buffer);
    BodyAlignedBoxByParticle emitter_b2(oil, makeShared<AlignedBoxShape>(xAxis, inletb2_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_b2(emitter_b2);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_b2(emitter_b2, inlet_buffer);
    BodyAlignedBoxByParticle emitter_b3(oil, makeShared<AlignedBoxShape>(xAxis, inletb3_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_b3(emitter_b3);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_b3(emitter_b3, inlet_buffer);
    BodyAlignedBoxByParticle emitter_b4(oil, makeShared<AlignedBoxShape>(xAxis, inletb4_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_b4(emitter_b4);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_b4(emitter_b4, inlet_buffer);
    BodyAlignedBoxByParticle emitter_b5(oil, makeShared<AlignedBoxShape>(xAxis, inletb5_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition_b5(emitter_b5);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_b5(emitter_b5, inlet_buffer);
    BodyAlignedBoxByCell outlet_disposer_a(oil, makeShared<AlignedBoxShape>(xAxis, outleta_transform, outlet_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> outlet_disposer_outflow_deletion_a(outlet_disposer_a);
    BodyAlignedBoxByCell outlet_disposer_b(oil, makeShared<AlignedBoxShape>(xAxis, outletb_transform, outlet_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> outlet_disposer_outflow_deletion_b(outlet_disposer_b);

    IsotropicDiffusion diffusion_oil("Phi", "Phi", k_oil);
    IsotropicDiffusion diffusion_winding("Phi", "Phi", k_winding);
    IsotropicDiffusion conductivity_oil_winding("Phi", "Phi", k_contact);
    ThermalRelaxationComplex thermal_relaxation_complex_oil(
        ConstructorArgs(oil_inner, &diffusion_oil),
        ConstructorArgs(oil_thermal_contact, &conductivity_oil_winding));
    ThermalRelaxationComplex thermal_relaxation_complex_winding(
        ConstructorArgs(winding_inner, &diffusion_winding),
        ConstructorArgs(winding_thermal_contact, &conductivity_oil_winding));
    SimpleDynamics<ThermoWindingInitialCondition> thermalwinding_condition(winding);
    SimpleDynamics<ThermofluidBodyInitialCondition> thermalfluid_initial_condition(oil);
    SimpleDynamics<ThermoWindingHeatSource> heat_source(winding);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(oil, "Indicator");
    body_states_recording.addToWrite<Real>(oil, "Phi");
    body_states_recording.addToWrite<Real>(winding, "Phi");
    // RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(oil, gravity);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_slot_phi("Phi", winding_slot_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_tooth_phi("Phi", winding_tooth_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_yoke_phi("Phi", winding_yoke_contact);
    //----------------------------------------------------------------------
    //	Building of multibody for rotor rotation.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the rigid shell box. */
    SolidBodyPartForSimbody rotor_multibody(rotor, makeShared<RotorFromMesh>("Rotor"));
    SimTK::Body::Rigid rigid_info(*rotor_multibody.body_part_mass_properties_);
    SimTK::MobilizedBody::Pin
        Rotor_Pin(matter.Ground(), SimTK::Transform(SimTKVec3(0)), rigid_info, SimTK::Transform(SimTKVec3()));
    /** Initial angle of rotation. */
    Rotor_Pin.setDefaultAngle(0.0);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    Rotor_Pin.setRate(state, Omega);
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.initialize(state);
    /** Coupling between SimBody and SPH.*/
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_rotor(rotor_multibody, MBsystem, Rotor_Pin, integ);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    hausing_normal_direction.exec();
    rotor_normal_direction.exec();
    winding_normal_direction.exec();
    cover_normal_direction.exec();
    thermalwinding_condition.exec();
    thermalfluid_initial_condition.exec();
    indicate_free_surface.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Time stepping control parameters.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 30.0;
    Real output_interval = 0.05;
    Real dt = 0.0; /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    // write_water_mechanical_energy.writeToFile(number_of_iterations);
    write_slot_phi.writeToFile(number_of_iterations);
    write_tooth_phi.writeToFile(number_of_iterations);
    write_yoke_phi.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                SimTK::State &state_for_update = integ.updAdvancedState();
                integ.stepBy(dt);
                constraint_rotor.exec();

                pressure_relaxation.exec(dt);
                inflow_condition_a1.exec();
                inflow_condition_a2.exec();
                inflow_condition_a3.exec();
                inflow_condition_a4.exec();
                inflow_condition_a5.exec();
                inflow_condition_b1.exec();
                inflow_condition_b2.exec();
                inflow_condition_b3.exec();
                inflow_condition_b4.exec();
                inflow_condition_b5.exec();
                density_relaxation.exec(dt);
                inflow_condition_a1.exec();
                inflow_condition_a2.exec();
                inflow_condition_a3.exec();
                inflow_condition_a4.exec();
                inflow_condition_a5.exec();
                inflow_condition_b1.exec();
                inflow_condition_b2.exec();
                inflow_condition_b3.exec();
                inflow_condition_b4.exec();
                inflow_condition_b5.exec();
                thermal_relaxation_complex_oil.exec(dt);
                heat_source.exec(dt); /** Implement of heat recources. */
                thermal_relaxation_complex_winding.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                TickCount current_time = TickCount::now();
                TimeInterval elapsed_time = current_time - t1 - interval;
                Real remaining_physical_time = end_time - GlobalStaticVariables::physical_time_;
                Real estimated_remaining_real_time = elapsed_time.seconds() *
                                                     (remaining_physical_time / GlobalStaticVariables::physical_time_);
                int hours = static_cast<int>(estimated_remaining_real_time) / 3600;
                int minutes = (static_cast<int>(estimated_remaining_real_time) % 3600) / 60;
                int seconds = static_cast<int>(estimated_remaining_real_time) % 60;
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << " Time = " << GlobalStaticVariables::physical_time_
                          << " Dt = " << Dt << " dt = " << dt
                          << " Remaining: " << hours << " h "
                          << minutes << " min " << seconds << " s\n";
            }
            number_of_iterations++;

            /** inflow emitter injection*/
            emitter_injection_a1.exec();
            emitter_injection_a2.exec();
            emitter_injection_a3.exec();
            emitter_injection_a4.exec();
            emitter_injection_a5.exec();
            emitter_injection_b1.exec();
            emitter_injection_b2.exec();
            emitter_injection_b3.exec();
            emitter_injection_b4.exec();
            emitter_injection_b5.exec();
            /** outflow delete*/
            outlet_disposer_outflow_deletion_a.exec();
            outlet_disposer_outflow_deletion_b.exec();
            /** Update cell linked list and configuration. */
            oil.updateCellLinkedListWithParticleSort(100);
            oil_body_complex.updateConfiguration();
            oil_thermal_complex.updateConfiguration();
            winding_thermal_complex.updateConfiguration();
            winding_slot_contact.updateConfiguration();
            winding_tooth_contact.updateConfiguration();
            winding_yoke_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        // write_water_mechanical_energy.writeToFile(number_of_iterations);
        indicate_free_surface.exec();
        body_states_recording.writeToFile();
        write_slot_phi.writeToFile(number_of_iterations);
        write_tooth_phi.writeToFile(number_of_iterations);
        write_yoke_phi.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    double total_seconds = tt.seconds();
    int hours = static_cast<int>(total_seconds / 3600);
    int remaining = static_cast<int>(total_seconds) % 3600;
    int minutes = remaining / 60;
    int seconds = remaining % 60;
    std::cout << "Total wall time for computation: "
              << hours << " hours "
              << minutes << " minutes "
              << seconds << " seconds." << std::endl;

    return 0;
}
