/**
 * @file 	filling_tank_vertical.cpp
 * @brief 	2D example to show that a tank is verticaL filled by emitter.
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;              /**< Tank length. */
Real DH = 5.366;              /**< Tank height. */
Real resolution_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4; /**< Extending width for wall boundary. */
Real LL = 0.125;           /**< Inflow region length. */
Real LH = 2 * BW;              /**< Inflows region height. */
Real inlet_height = DH - 0.5 * LH; /**< Inflow location height */
Vec2d inlet_halfsize = Vec2d(0.5 * LH, 0.5 * LL);
Vec2d inlet_translation = Vec2d(0.5 * DL, DH);
Transform inlet_transform(Rotation2d( - (Pi / 2)), inlet_translation);
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
// Create a rotor.
Real DS = 2; /**< Diameter of transmission shaft on rotor. */
Real RS = 0.5 * DS; /**< Radius of transmission shaft on rotor. */
Real reference_circle_radius = 0.25 * RS; /**< Radius of reference circlev on rotor. */
Vecd center(DL / 2, DH / 2); /**< Location of transmission shaft on rotor. */
Vec2d reference_circle_initial_center = center + Vec2d(0.0, 0.5 * RS);
int resolution_circle = 100; /**<Approximate the circle as the number of sides of the polygon. */
Real rotor_angular_velocity = 5.0; /**<Angular velocity. */
Real rotor_current_angle = 0.0;    /**<Angle of rotor. */
Transform rotor_rotation_transform(Rotation2d(rotor_current_angle), center);
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
//	Rotated rotor
//----------------------------------------------------------------------
std::vector<Vecd> CreateRotatingCircleShape(Vecd center, Real radius, int resolution_circle, Real current_angle)
{
    std::vector<Vecd> circle_shape;

    for (int i = 0; i <= resolution_circle; ++i)
    {
        // Calculate the angle for each point on the circle
        Real theta = 2.0 * Pi * i / resolution_circle;

        // Initial position on the circle (without rotation)
        Real x = center[0] + radius * cos(theta);
        Real y = center[1] + radius * sin(theta);

        // Apply rotation to the point (using the current_angle)
        Real rotated_x = center[0] + (x - center[0]) * cos(current_angle) - (y - center[1]) * sin(current_angle);
        Real rotated_y = center[1] + (x - center[0]) * sin(current_angle) + (y - center[1]) * cos(current_angle);

        circle_shape.push_back(Vecd(rotated_x, rotated_y));
    }

    return circle_shape;
}
Vec2d CreateRotatingReferenceCircleCenter(Vecd initial_center, Vecd center, Real current_angle)
{
    Real rotated_x = center[0] + (initial_center[0] - center[0]) * cos(current_angle) - (initial_center[1] - center[1]) * sin(current_angle);
    Real rotated_y = center[1] + (initial_center[0] - center[0]) * sin(current_angle) + (initial_center[1] - center[1]) * cos(current_angle);
    return Vec2d(rotated_x, rotated_y);
}
//----------------------------------------------------------------------
//	Case-dependent wall boundary
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
    std::vector<Vecd> rotated_circle_shape = CreateRotatingCircleShape(center, RS, resolution_circle, rotor_current_angle);
    Vec2d rotated_reference_circle_center = CreateRotatingReferenceCircleCenter(reference_circle_initial_center, center, rotor_current_angle);
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(CreateOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(CreateInnerWallShape(), ShapeBooleanOps::sub);
        multi_polygon_.addABox(inlet_transform, inlet_halfsize, ShapeBooleanOps::sub);
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
        return Vec2d(0.5, 0.0);
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
    water_body.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(350.0);
    water_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    auto wall_boundary = makeShared<WallBoundary>("Wall");
    SolidBody wall(sph_system, wall_boundary);
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();
    
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
    ContactRelation water_body_contact(water_body, {&wall});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_body});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_body_complex(water_body_inner, water_body_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(water_body_inner, water_body_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_body_inner, water_body_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_body_inner, water_body_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_body_inner, water_body_contact);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_body, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);

    BodyAlignedBoxByParticle emitter(water_body, makeShared<AlignedBoxShape>(xAxis, inlet_transform , inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, inlet_buffer);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(water_body, "Indicator");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_body, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    indicate_free_surface.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Time stepping control parameters.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 30.0;
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

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);
                inflow_condition.exec();
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
            fluid_observer_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_water_mechanical_energy.writeToFile(number_of_iterations);
        indicate_free_surface.exec();
        body_states_recording.writeToFile();
        write_recorded_water_pressure.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    write_water_mechanical_energy.testResult();
    write_recorded_water_pressure.testResult();

    return 0;
}
