#pragma once
#include <core/Scene.hpp>
#include <pkg/common/Recorder.hpp>

namespace yade { // Cannot have #include directive inside.

class ForceRecorder : public Recorder {
public:
	void action() override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(ForceRecorder,Recorder,"Engine saves the resultant force affecting to bodies, listed in `ids`. For instance, can be useful for defining the forces, which affects to _buldozer_ during its work.",
		((std::vector<int>,ids,,,"List of bodies whose state will be measured"))
		((Vector3r,totalForce,Vector3r::Zero(),,"Resultant force, returning by the function."))
		,
		initRun=true;
		,
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ForceRecorder);

class TorqueRecorder : public Recorder {
public:
	void action() override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(TorqueRecorder,Recorder,"Engine saves the total torque according to the given axis and ZeroPoint, the force is taken from bodies, listed in `ids`  For instance, can be useful for defining the torque, which affects on ball mill during its work.",
		((std::vector<int>,ids,,,"List of bodies whose state will be measured"))
		((Vector3r,rotationAxis,Vector3r::UnitX(),,"Rotation axis"))
		((Vector3r,zeroPoint,Vector3r::Zero(),,"Point of rotation center"))
		((Real,totalTorque,0,,"Resultant torque, returning by the function.")),
		initRun=true;
	);
	// clang-format on
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(TorqueRecorder);

} // namespace yade
