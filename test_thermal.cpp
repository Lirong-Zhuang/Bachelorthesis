/**
 * @file 	filling_tank_vertical.cpp
 * @brief 	2D example to show that a tank is verticaL filled by emitter.
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                   /**< Tank length. */
Real DH = 5.366;                   /**< Tank height. */
Real resolution_ref = 0.025;       /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;      /**< Extending width for wall boundary. */
Real LL = 0.125;                   /**< Inflow region length. */
Real LH = 2 * BW;                  /**< Inflows region height. */
Real inlet_height = DH - 0.5 * LH; /**< Inflow location height */
Vec2d inlet_halfsize = Vec2d(0.5 * LH, 0.5 * LL);
Vec2d inlet_translation = Vec2d(0.5 * DL, DH);
Transform inlet_transform(Rotation2d(-(Pi / 2)), inlet_translation);
Real inlet_height2 = 1.0;   /**< Inflow location height */
Real inlet_distance2 = -BW; /**< Inflow location distance */
Vec2d inlet_halfsize2 = Vec2d(0.5 * LH, 0.5 * LL);
Vec2d inlet_translation2 = Vec2d(inlet_distance2, inlet_height2) + inlet_halfsize2;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
Real rho0_f = 1.0;                                      /**< Reference density of fluid. */
Real gravity_g = 1.0;                                   /**< Gravity force of fluid. */
Real U_f = 2.0 * sqrt(gravity_g * (inlet_height + LH)); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                  /**< Reference sound speed. */
Real diffusion_coeff = 1.0e-3;
Real diffusion_coeff_rotor = 1.0e-3;
Real K = (diffusion_coeff * diffusion_coeff_rotor) / (diffusion_coeff + diffusion_coeff_rotor);
    Real Re = 100.0;                    /**< Reynolds number100. */
Real mu_f = rho0_f * U_f * LL / Re; /**< Dynamics viscosity. */
// Create a rotor.
Real DS = 2;                                           /**< Diameter of transmission shaft on rotor. */
Real RS = 0.5 * DS;                                    /**< Radius of transmission shaft on rotor. */
Real reference_circle_radius = 0.25 * RS;              /**< Radius of reference circlev on rotor. */
Vecd center(DL / 2, DH / 2);                           /**< Location of transmission shaft on rotor. */
int resolution_circle = 50;                            /**<Approximate the circle as the number of sides of the polygon. */
Real rotor_rotation_velocity = 600;                    /**<Angular velocity rpm. */
Real Omega = -(rotor_rotation_velocity * 2 * Pi / 60); /**<Angle of rotor. */
//----------------------------------------------------------------------
//	Global parameters on the initial condition
//----------------------------------------------------------------------
Real phi_rotor = 100.0;
Real phi_fluid_initial = 0.0;
//----------------------------------------------------------------------
//	Geometries
//----------------------------------------------------------------------
/** create a outer wall polygon. */
std::vector<Vecd> CreateOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-BW, -BW));

    return outer_wall_shape;
}
/** create a inner wall polygon. */
std::vector<Vecd> CreateInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, DH));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, 0.0));

    return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Case-dependent wall boundary
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(CreateOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(CreateInnerWallShape(), ShapeBooleanOps::sub);
        multi_polygon_.addABox(inlet_transform, inlet_halfsize, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Application dependent initial condition of rotor.
//----------------------------------------------------------------------
class ThermoRotorInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermoRotorInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};
    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = phi_rotor;
    };

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
//----------------------------------------------------------------------
//	Case-dependent rotor boundary
//----------------------------------------------------------------------
class RotorBoundary : public MultiPolygonShape
{
  public:
    explicit RotorBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, RS, resolution_circle, ShapeBooleanOps::add);
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
        return Vec2d(3.0, 0.0);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Build up a SPHSystem */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    TransformShape<GeometricShapeBox> water_inlet_shape(inlet_transform, inlet_halfsize);
    FluidBody water_body(sph_system, water_inlet_shape, "WaterBody");
    water_body.sph_adaptation_->resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    water_body.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(350.0);
    water_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    SolidBody rotor(sph_system, makeShared<RotorBoundary>("Rotor"));
    rotor.defineMaterial<Solid>();
    rotor.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_body);
    InnerRelation rotor_inner(rotor);
    ContactRelation water_body_contact(water_body, {&wall, &rotor});
    ContactRelation water_rotor_contact(water_body, {&rotor});
    ContactRelation rotor_contact(rotor, {&water_body});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_body});
    ContactRelation rotor_observer_contact(fluid_observer, {&rotor});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_body_complex(water_body_inner, water_body_contact);
    ComplexRelation water_rotor_complex(water_body_inner, water_rotor_contact);
    ComplexRelation rotor_complex(rotor_inner, rotor_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    SimpleDynamics<NormalDirectionFromBodyShape> rotor_normal_direction(rotor);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(water_body_inner, water_body_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_body_inner, water_body_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_body_inner, water_body_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_body_inner, water_body_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_body_inner, water_body_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(water_body_inner, water_body_contact);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_body, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);

    BodyAlignedBoxByParticle emitter(water_body, makeShared<AlignedBoxShape>(xAxis, inlet_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, inlet_buffer);

    IsotropicDiffusion diffusion("Phi", "Phi", diffusion_coeff);
    IsotropicDiffusion contact("Phi", "Phi", K);
    ThermalRelaxationComplex thermal_relaxation_complex(
        ConstructorArgs(water_body_inner, &diffusion),
        ConstructorArgs(water_rotor_contact, &contact));
    IsotropicDiffusion diffusion_rotor("Phi", "Phi", diffusion_coeff_rotor);
    ThermalRelaxationComplex thermal_relaxation_complex_rotor(
        ConstructorArgs(rotor_inner, &diffusion_rotor),
        ConstructorArgs(rotor_contact, &contact));
    SimpleDynamics<ThermoRotorInitialCondition> thermorotor_condition(rotor);
    SimpleDynamics<ThermofluidBodyInitialCondition> thermofluid_initial_condition(water_body);

    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_body_inner);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(water_body, "Indicator");
    body_states_recording.addToWrite<Real>(water_body, "Phi");
    body_states_recording.addToWrite<Real>(rotor, "Phi");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_body, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_fluid_phi("Phi", fluid_observer_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_rotor_phi("Phi", rotor_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    rotor_normal_direction.exec();
    thermorotor_condition.exec();
    thermofluid_initial_condition.exec();
    indicate_free_surface.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Time stepping control parameters.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 50.0;
    Real output_interval = 0.1;
    Real dt = 0.0; /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
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
            transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);
                inflow_condition.exec();
                thermal_relaxation_complex.exec(dt);
                thermal_relaxation_complex_rotor.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** inflow emitter injection*/
            emitter_injection.exec();
            /** Update cell linked list and configuration. */

            water_body.updateCellLinkedListWithParticleSort(100);
            water_body_complex.updateConfiguration();
            water_rotor_complex.updateConfiguration();
            rotor_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            rotor_observer_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        compute_vorticity.exec();
        write_water_mechanical_energy.writeToFile(number_of_iterations);
        indicate_free_surface.exec();
        body_states_recording.writeToFile();
        write_recorded_water_pressure.writeToFile(number_of_iterations);
        write_fluid_phi.writeToFile(number_of_iterations);
        write_rotor_phi.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
