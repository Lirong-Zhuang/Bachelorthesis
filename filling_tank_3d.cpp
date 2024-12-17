/**
 * @file 	filling_tank_3D.cpp
 * @brief 	3D model filling tank.
 * @author 	Lirong Zhuang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry, material parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;                   /**< Tank length. */
Real DH = 5.0;                   /**< Tank height. */
Real DW = 5.0;                   /**< Tank width. */
// Axis
Vec3d axis_x(1.0, 0.0, 0.0);
Vec3d axis_y(0.0, 1.0, 0.0);
Vec3d axis_z(0.0, 0.0, 1.0);
Real resolution_ref = 0.025;       /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;      /**< Extending width for wall boundary. */
Real LL = 0.25;                   /**< Inflow region diameter. */
Real LH = 2 * BW;                  /**< Inflows region height. */
Real OL = 0.25;                   /**< Inflow region diameter. */
Real OH = BW;                  /**< Inflows region height. */
Real inlet_height = DH - 0.5 * LH; /**< Inflow location height */
Vecd inlet_halfsize = Vecd(0.5 * LH, 0.5 * LL, 0.5 * LL);
Vecd inlet_translation = Vecd(0.0, 0.5 * DH + BW - 0.5 * LH, 0.0);
Transform inlet_transform(Rotation3d((-(Pi / 2)), axis_z), inlet_translation);
Vecd outlet_halfsize = Vecd(0.5 * OL, 0.5 * OH, 0.5 * OL);
Vecd outlet_translation = Vecd(0.5, -0.5 * DH - BW + 0.5 * OH, 0.0);
Transform outlet_transform(Rotation3d((-(Pi / 2)), axis_z), outlet_translation);
Real inlet_height2 = 1.0;   /**< Inflow location height */
Real inlet_distance2 = -BW; /**< Inflow location distance */
Vecd inlet_halfsize2 = Vecd(0.5 * LH, 0.5 * LL, 0.5 * LL);
Vecd inlet_translation2 = Vecd(-(0.5*DL+BW-0.5*LH), 0.0, 0.0);
Real rho0_f = 1.0;                                      /**< Reference density of fluid. */
Real gravity_g = 1.0;                                   /**< Gravity force of fluid. */
Real U_f = 2.0 * sqrt(gravity_g * (inlet_height + LH)); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                  /**< Reference sound speed. */
Real Re = 100.0;     /**< Reynolds number100. */
Real mu_f = 0.00975; // rho0_f * U_f * LL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-0.5 * DL - 2 * BW, -0.5 * DH - 2 * BW, -0.5 * DW - 2 * BW);
Vec3d domain_upper_bound(0.5 * DL + 2 * BW, 0.5 * DH + 2 * BW, 0.5 * DW + 2 * BW);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
// dynamics informations of rotor
Real rotor_rotation_velocity = 900;                    /**<Angular velocity rpm. */
Real Omega = -(rotor_rotation_velocity * 2 * Pi / 60); /**<Angle of rotor. */
// thermal parameters
Real phi_winding = 110.0;                                       /**< Temperature of winding at begin. */
Real phi_fluid_initial = 75.0;                                  /**< Temperature of oil at begin. */
Real k_oil = 7.62e-8;                                           /**< Diffusion coefficient of oil 2.0e-7. */
Real k_winding = 1.14e-4;                                       /**< Diffusion coefficient of winding 1.1e-6. */
Real k_contact = (2 * k_oil * k_winding) / (k_oil + k_winding); /**< Thermal conductivity between winding and oil. */
Real dq = 1.5;                                                  /**< Heating efficient of internal heat source [°C/s]. */
//----------------------------------------------------------------------
//	Case-dependent Tank
//----------------------------------------------------------------------
class Wallboundary : public ComplexShape
{
  public:
    explicit Wallboundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Vecd translation_tank(0.0, 0.0, 0.0);
        add<TriangleMeshShapeBrick>(halfsize_outer, resolution_ref, translation_tank);
        subtract<TriangleMeshShapeBrick>(halfsize_inner, resolution_ref, translation_tank);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), 0.5 * LL, 0.5 * LH, resolution_ref, inlet_translation);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0, 0.0), 0.5 * LL, 0.5 * LH, resolution_ref, inlet_translation2);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), 0.5 * OL, 0.5 * OH, resolution_ref, outlet_translation);
    }
};
    //----------------------------------------------------------------------
    //	Case-dependent Fluid
    //----------------------------------------------------------------------
    class Fluidboundary : public ComplexShape
    {
      public:
        explicit Fluidboundary(const std::string &shape_name) : ComplexShape(shape_name)
        {
            add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 1.0, 0), 0.5 * LL, 0.5 * LH, resolution_ref, inlet_translation);
            // add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0, 0, 1.0), 0.5 * LL, 0.5 * LH, resolution_ref, inlet_translation2);
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
            return Vec3d(1.0, 0.0, 0.0);
        }
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
        /** Set the starting time. */
        GlobalStaticVariables::physical_time_ = 0.0;
        //----------------------------------------------------------------------
        //	Creating body, materials and particles.
        //----------------------------------------------------------------------
        SolidBody tank(sph_system, makeShared<Wallboundary>("Tank"));
        tank.defineMaterial<Solid>();
        tank.generateParticles<BaseParticles, Lattice>();

        FluidBody oil(sph_system, makeShared<Fluidboundary>("Oil"));
        oil.sph_adaptation_->resetKernel<KernelTabulated<KernelWendlandC2>>(20);
        oil.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
        ParticleBuffer<ReserveSizeFactor> inlet_buffer(3000.0);
        oil.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);
        //----------------------------------------------------------------------
        //	Define body relation map.
        //	The contact map gives the topological connections between the bodies.
        //	Basically the the range of bodies to build neighbor particle lists.
        //  Generally, we first define all the inner relations, then the contact relations.
        //  At last, we define the complex relaxations by combining previous defined
        //  inner and contact relations.
        //----------------------------------------------------------------------
        InnerRelation oil_inner(oil);
        ContactRelation oil_solid_contact(oil, {&tank});
        //----------------------------------------------------------------------
        // Combined relations built from basic relations
        //----------------------------------------------------------------------
        ComplexRelation oil_body_complex(oil_inner, oil_solid_contact);
        //----------------------------------------------------------------------
        //	Define all numerical methods which are used in this case.
        //----------------------------------------------------------------------
        SimpleDynamics<NormalDirectionFromBodyShape> tank_normal_direction(tank);
        InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> indicate_free_surface(oil_inner, oil_solid_contact);

        Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(oil_inner, oil_solid_contact);
        Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(oil_inner, oil_solid_contact);
        InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(oil_inner, oil_solid_contact);

        Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
        SimpleDynamics<GravityForce> constant_gravity(oil, gravity);
        ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(oil, U_f);
        ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(oil);
        InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(oil_inner, oil_solid_contact);

        BodyAlignedBoxByParticle emitter(oil, makeShared<AlignedBoxShape>(xAxis, inlet_transform, inlet_halfsize));
        SimpleDynamics<InletInflowCondition> inflow_condition(emitter);
        SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection(emitter, inlet_buffer);
        BodyAlignedBoxByCell outlet_disposer(oil, makeShared<AlignedBoxShape>(xAxis, outlet_transform, outlet_halfsize));
        SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> outlet_disposer_outflow_deletion(outlet_disposer);
        //----------------------------------------------------------------------
        //	File output and regression check.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp body_states_recording(sph_system);
        body_states_recording.addToWrite<int>(oil, "Indicator");
        RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(oil, gravity);
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        tank_normal_direction.exec();
        indicate_free_surface.exec();
        constant_gravity.exec();
        //----------------------------------------------------------------------
        //	Time stepping control parameters.
        //----------------------------------------------------------------------
        size_t number_of_iterations = sph_system.RestartStep();
        int screen_output_interval = 100;
        Real end_time = 20.0;
        Real output_interval = 0.05;
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
                emitter_injection.exec();
                /** outflow delete*/
                outlet_disposer_outflow_deletion.exec();
                /** Update cell linked list and configuration. */
                oil.updateCellLinkedListWithParticleSort(100);
                oil_body_complex.updateConfiguration();
            }

            TickCount t2 = TickCount::now();
            write_water_mechanical_energy.writeToFile(number_of_iterations);
            indicate_free_surface.exec();
            body_states_recording.writeToFile();
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
