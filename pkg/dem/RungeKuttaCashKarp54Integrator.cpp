#ifdef YADE_ODEINT
#include <core/Scene.hpp>
#include <pkg/dem/RungeKuttaCashKarp54Integrator.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((RungeKuttaCashKarp54Integrator));

shared_ptr<Integrator> RungeKuttaCashKarp54Integrator_ctor_list(const boost::python::list& slaves)
{
	shared_ptr<Integrator> instance(new RungeKuttaCashKarp54Integrator);
	instance->slaves_set(slaves);
	return instance;
}

void RungeKuttaCashKarp54Integrator::action()
{
	Real dt = scene->dt;

	Real time = scene->time;

	// declaration of ‘rungekuttaerrorcontroller’ shadows a member of ‘yade::RungeKuttaCashKarp54Integrator’ [-Werror=shadow]
	error_checker_type rungekuttaerrorcontroller2 = error_checker_type(abs_err, rel_err, a_x, a_dxdt);

	// declaration of ‘rungekuttastepper’ shadows a member of ‘yade::RungeKuttaCashKarp54Integrator’ [-Werror=shadow]
	controlled_stepper_type rungekuttastepper2 = controlled_stepper_type(rungekuttaerrorcontroller2);

	stateVector currentstates = getCurrentStates();

	resetstate.reserve(currentstates.size());

	copy(currentstates.begin(), currentstates.end(), back_inserter(resetstate)); //copy current state to resetstate

	this->timeresetvalue = time; //set reset time to the time just before the integration

	/*Try an adaptive integration*/

	integrationsteps += integrate_adaptive(
	        rungekuttastepper2, make_ode_wrapper(*((Integrator*)this), &Integrator::system), currentstates, time, time + dt, stepsize, observer(this));

	scene->time = scene->time - dt; //Scene move next time step function already increments the time so we have to decrement it just before it.
}

} // namespace yade

#endif
