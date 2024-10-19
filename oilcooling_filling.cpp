/**
 * @file 	oilcooling_em.cpp
 * @brief 	2D example to show the filling process of oil cooling system in a running electric motor.
 * @details This is the feasibility verification for the 3D model based on the experimental study of the team
 * 			of Tanguy Davin.Coding ist based on the programm filling_tank.cpp of Professor Xiangyu Hu at TUM.
 * @author 	Lirong Zhuang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DM = 200 / 10; /**< Diameter of motor hausing. */
Real DR = 120 / 10; /**< Diameter of rotor. */
Real DS = 40 / 10;  /**< Diameter of transmission shaft on rotor. */
Real WL = 30 / 10;  /**< Winding length. */
Real WH = 25 / 10;  /**< Winding height. */
Real AG = 1 / 10;   /**< Air-gap. */
Real LL = 2.8 / 10; /**< Inflow region length. */

Real PI = 3.1415926;
Real RM = 0.5 * DM;                      /**< Radius of motor hausing. */
Real RR = 0.5 * DR;                      /**< Radius of rotor. */
Real RS = 0.5 * DS;                      /**< Radius of transmission shaft on rotor. */
Real Wnum = 12;                          /**< Winding Number. */
Real WD = RR + AG + 0.5 * WH;            /**< Distance from center point to the center of a winding. */
Real ZW = RR + AG;                       /**< Distance from center point to the site of a winding. */
int resolution_circle = 200;             /**<Approximate the circle as the number of sides of the polygon. */
Real resolution_ref = 0.025;             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;            /**< Extending width for wall boundary. */
Real LH = 2.0 * BW;                      /**< Inflows region height. */
Real Lnum = 5;                           /**< Inflows number. */
Real Lstart = (Wnum - 2 * Lnum + 2) / 4; /**< Start location of inlets. */
Real Lend = (Wnum + 2 * Lnum - 2) / 4;   /**< End location of inlets. */
Real inlet_height = RM - BW;             /**< Inflow location height */
Real inlet_distance = -(0.5 * LL);          /**< Inflow location distance */
Vec2d inlet_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d inlet_translation = Vec2d(inlet_distance, inlet_height) + inlet_halfsize;
Transform inlet_transform(Rotation2d(-(Pi / 2)), inlet_translation);
BoundingBox system_domain_bounds(Vec2d(-BW - RM, -BW - RM), Vec2d(RM + BW, RM + BW));
Vecd center(0.0, 0.0);
// observer location
StdVec<Vecd> observation_location = {Vecd(0, 0)};
Real rho0_f = 0.9;                                      /**< Reference density of fluid. */
Real gravity_g = 1.0;                                   /**< Gravity force of fluid. */
Real U_f = 2.0 * sqrt(gravity_g * (inlet_height + LH)); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                  /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometries
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Case-dependent wall boundary
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, (RM + BW) , resolution_circle, ShapeBooleanOps::add);         /**< Outer wall of motor hausing. */
        multi_polygon_.addACircle(center, RM , resolution_circle, ShapeBooleanOps::sub);                /**< Inner wall of motor hausing. */
        multi_polygon_.addABox(inlet_transform, inlet_halfsize, ShapeBooleanOps::sub);     /**< Inlets. */
        multi_polygon_.addACircle(center, RS , resolution_circle, ShapeBooleanOps::add);                /**< Shaft. */
        /** Add the windings. */
        Real angle_increment = 2 * PI / Wnum;
        for (int i = 0; i < Wnum; ++i)
        {
            Real theta = i * angle_increment;
            Real center_x = WD * cos(theta);
            Real center_y = WD * sin(theta);
            Vec2d winding_translation(center_x, center_y);
            Vec2d winding_halfsize(WL / 2 , WH / 2);
            Transform winding_transform(Rotation2d(theta - (PI / 2)), winding_translation);
            multi_polygon_.addABox(winding_transform, winding_halfsize, ShapeBooleanOps::add);
        }
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
        return Vec2d(1.0, 0.0);
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
    TransformShape<GeometricShapeBox> oil_inlet_shape(inlet_transform, inlet_halfsize);
    FluidBody oil_body(sph_system, oil_inlet_shape, "OilBody");
    oil_body.sph_adaptation_->resetKernel<KernelTabulated<KernelWendlandC2>>(20);
    oil_body.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(3500.0);
    oil_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    SolidBody wall(sph_system, makeShared<WallBoundary>("Wall"));
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
    InnerRelation oil_body_inner(oil_body);
    ContactRelation oil_body_contact(oil_body, {&wall});
    ContactRelation fluid_observer_contact(fluid_observer, {&oil_body});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation oil_body_complex(oil_body_inner, oil_body_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(oil_body_inner, oil_body_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(oil_body_inner, oil_body_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(oil_body_inner, oil_body_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(oil_body_inner, oil_body_contact);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(oil_body, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(oil_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(oil_body);

    BodyAlignedBoxByParticle emitter(oil_body, makeShared<AlignedBoxShape>(xAxis, inlet_transform, inlet_halfsize));
    SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, inlet_buffer);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(oil_body, "Indicator");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(oil_body, gravity);
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

            oil_body.updateCellLinkedListWithParticleSort(100);
            oil_body_complex.updateConfiguration();
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
