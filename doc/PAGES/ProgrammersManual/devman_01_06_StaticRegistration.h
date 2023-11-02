/**\page DevManStaticRegistration Static Registration

\tableofcontents

This page explains the details of the two utility classes `ParameterBlock` and
`InputParameters`.

\section DevManParametersSec1 1 The general philosophy of Static Registration

First, what is <B>Static Registration</B>? When a c++ program is executed there
are a number of phases that occur in order. Once such phase, before entering the
`main` function, is the phase in which all static variables are initialized. We
shall call this phase the static initialization phase.
In the most basic sense static registration is the process whereby we create a
registry of object-constructors during this phase. This registry can then used,
during the dynamic phases of the program, to construct objects.
A simple example is the code below:
\code
#include <iostream>
#include "DoSomethingFunctions"

static int a = DoSomething();

int main()
{
  std::cout << "The value of a is: " << a << std::endl;
}
\endcode
where the variable `int a` is initialized, before the `main` function starts,
using a function call, which happens during the static initialization phase of
the program.

\subsection DevManParametersSec1_1 Rigging static initialization for object factories or object makers

We can use the rather odd dynamics of singleton objects to create a static
registry of function creation calls during static initialization. In the code
below we define an `ObjectMaker` class that is a singleton. A singleton is an
object of which only a single instance can be created, which we do by declaring
a static instance, a bunch of static members methods, a private constructor and
deleted copy- and move constructors.

The `ObjectMaker` maintains a map of string keys to object construction calls.
The construction calls always return a pointer to a base class `Object` but can
construct any child object by using the shown template arguments. In order for
us to achieve basic static registration we need two templated-methods, these
are `CallObjectConstructor` and `RegisterObject` in the code below.
The first method defines a function-address that we can put into our registry.
The second method is the actual method we call to make a registry entry. When we
then want to make a registered object we simply call the `MakeRegisteredObject`
method.

\code
#include <iostream>
#include <map>
#include <functional>

class Object
{
public:
  int member_;
  Object() : member_(99) {}
  Object(int member) : member_(member) {}
};

class ChildObject : public Object
{
public:
  ChildObject() : Object(999) {}
};

//===============================================
typedef std::function<Object*()> CreationFunction;
class ObjectMaker // Singleton
{
public:
  std::map<std::string, CreationFunction> registry_;

public:
  static ObjectMaker& GetInstance() noexcept
  {
    static ObjectMaker singleton;
    return singleton;
  }
  template <typename T>
  static Object* CallObjectConstructor() { return new T(); }

  template <typename T>
  static int RegisterObject(const std::string& register_name)
  {
    auto& instance = GetInstance();
    instance.registry_.insert(
      std::make_pair(register_name, &CallObjectConstructor<T>));
    return 0;
  }

  static Object* MakeRegisteredObject(const std::string& obj_name)
  {
    auto& instance = GetInstance();

    return instance.registry_.at(obj_name)();
  }
private:
  ObjectMaker() = default;
  ObjectMaker(const ObjectMaker&) = delete;
  ObjectMaker(const ObjectMaker&&) = delete;
  ObjectMaker& operator=(const ObjectMaker&) = delete;
};

//===============================================
// This is where we statically register an object
static int abc12345 = ObjectMaker::RegisterObject<Object>("Parent");
static int abc12346 = ObjectMaker::RegisterObject<ChildObject>("Child");


int main()
{
  Object objectA(98); // Conventional object creation
  auto objectB = ObjectMaker::MakeRegisteredObject("Parent");
  auto objectC = ObjectMaker::MakeRegisteredObject("Child");

  std::cout << "Object A.member_ = " << objectA.member_ << "\n";
  std::cout << "Object B.member_ = " << objectB->member_ << "\n";
  std::cout << "Object C.member_ = " << objectC->member_ << "\n";

  return 0;
}
\endcode
Output:
\verbatim
Object A.member_ = 98
Object B.member_ = 99
Object C.member_ = 999
\endverbatim

Notice that we abuse the static initialization phase right above the `main`
function. We call the static singleton `ObjectMaker` and `RegisterObject` to
create registry entries for a `Parent` and a `Child`. In the main function we
then construct some of these statically registered objects.

\section DevManParametersSec2 2 Implementation in ChiTech: ChiObject and ChiObjectFactory
In ChiTech we have to, of course, enhance the functionality by a lot, to this
end we have the base object, `ChiObject`, and the factory, `ChiObjectFactory`.
Also, instead of using raw pointers, we use shared pointers for all
objects that can be statically registered.

\image html "DeveloperManual/StaticRegistration/StaticRegistration.png" "The factory pattern" width=900px

We developed a separate system to
control object construction parameters (more on this later). Additionally,
we do not want developers to come up with unique dummy variable names for the
static registration, i.e.,
`static int abc12345` or `static int abc12346`, just to make use of the
`ChiObjectFactory`, therefore, we use a macro. That works as follows.

In order to register an object within the `ChiObjectFactory` we created a macro
called `RegisterChiObject(namespace_name, obj_name)`. The arguments to this
macro, `namespace_name` and `obj_name`, then get used to construct the following
\code
RegisterChiObject(namespace_name, ObjectName);
//
// becomes:
//
static char unique_var_name_object_ObjectName_0 =
  ChiObjectFactory::AddObjectToRegistry<ObjectName, ChiObject>("namespace_name",
                                                             "ObjectName")
\endcode
which allows us to abuse static initialization to register the object with
text-name `"namespace_name::ObjectName"`. Note that nested namespaces are
allowed, e.g., `"space1::space1b::space1bc"`, is allowed.

We also implemented a number of "tricks". For example we created a
behind-the-scenes connection to the lua console so that registered objects can
be created from the lua input system. We also created a "RegistryDump" where
each registered object can dump its input arguments specification so we can use
it to build documentation (lessening the burden on the developer).

\section DevManParametersSec3 3 The InputParameters system
All statically registered objects have a common input parameters interface
implemented in the general class `chi::InputParameters`. This allows
us to fully support static registration if an object has the following basic
structure
\code
#include "ChiObject/chi_object.h"

namespace zorba
{

class CoolObject : public ChiObject
{
  const int a_;
  const double b_;
  //Misc. members
public:
  static chi::InputParameters GetInputParameters();
  CoolObject(const chi::InputParameters& params);

  //Misc. methods
};

}
\endcode
Here we have the static member method `GetInputParameters` which returns a
`chi::InputParameters` object used to assign parameters from user input,
and the constructor `CoolObject(const chi::InputParameters& params);`
which uses such an object.

In the definition (`.cc`) file we then have
\code
#include "CoolObject.h"

#include "ChiObject/object_maker.h"

namespace zorba
{

RegisterChiObject(chi_unit_tests, CoolObject);

chi::InputParameters CoolObject::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("A Cool object to play with");
  params.SetDocGroup("doc_Coolstuff"); // Documentation link

  params.AddRequiredParameter<int>("a", "Coolness factor a");
  params.AddOptionalParameter("b", 1.0, "Coolness factor b");

  return params;
}

CoolObject::CoolObject(const chi::InputParameters& params) :
  ChiObject(params),
  a_(params.GetParamValue<int>("a")),
  b_(params.GetParamValue<double>("b"))
{
}

}//namespace zorba
\endcode

Creating this object from lua is then as easy as
\code
params = {a = 98, b = 99.0}
cool_obj = zorba.CoolObject.Create(params)
\endcode

\subsection DevManParametersSecFAQ1 FAQ1 Documentation-strings (Doc-strings)
All parameters must have a documentation string (short name "doc string").
For example:
\code
params.AddRequiredParameter<int>("a", "Coolness factor a");
params.AddOptionalParameter("b", 1.0, "Coolness factor b");
\endcode
here "Coolness factor a" and "Coolness factor b" are the doc strings for
parameters "a" and "b" respectively.

You can also add documentation for the entire object. When the documentation
is auto-generated this description will appear at the top.
\code
params.SetGeneralDescription("General test object");
\endcode

By default the object's documentation will be posted to a group that is
reference-able via doxygen as `<namespace_name>__<ObjectName>` (note the double
underscore is a replacement for the namespace scope clarifiers `::`). You can
link the documentation to other groups by using the `SetDocGroup()`
method. This functionality allows you to customize your documentation according to your
needs.

<B>Note:</B> The block of text supplied to `SetGeneralDescription` is pasted
into a documentation page that will be processed by doxygen so all doxygen
commands will work here (including formulas, code blocks, and links).

\subsection DevManParametersSecFAQ2 FAQ2 Optional Parameters
Optional parameters must always have a default value. However, you do not
necessarily have to use it. In cases where you use a parameter only when
supplied,
the `InputParameters` object has the interface method `ParametersAtAssignment`
which returns the parameters used at assignment. It is used as follows:
\code
CoolObject::CoolObject(const chi::InputParameters& params) :
  ChiObject(params),
  a_(params.GetParamValue<int>("a")),
  b_(0.0)
{
  const auto& user_params = params.ParametersAtAssignment();
  if (user_params.Has("b"))
    b_ = user_params.GetParamValue<double>("b");
}
\endcode

\subsection DevManParametersSecFAQ3 FAQ3 Inheriting parameters from parent objects
Input parameters can be inherited from parents in the `GetInputParameters`
method call. For example
\code
chi::InputParameters LBSSolver::GetInputParameters()
{
  chi::InputParameters params =
    chi_physics::Solver::GetInputParameters();
  params += chi::SecondParent::GetInputParameters();
    //
    //etc.
\endcode
Gets all the parameters from the parent object `chi_physics::Solver` and then
`chi::SecondParent`, allowing us to change or add to these parameters.

In fact, a pretty slick move here is that we are not bound by classical
inheritance. The input parameters of any statically registered object can be
snatched and used by simply calling its static `GetInputParameters` method.

\subsection DevManParametersSecFAQ4 FAQ4 Changing inherited options to optional/required
If a required parameter is to be changed to optional we can do this using the
interface `ChangeExistingParamToOptional`, e.g.,
\code
chi::InputParameters LBSSolver::GetInputParameters()
{
  chi::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "LBSDatablock");
\endcode

Vice versa there is also `ChangeExistingParamToRequired<type>("name")`.

\subsection DevManParametersSecFAQ5 FAQ5 Deprecated parameters
Parameters in an `InputParameters` block can be marked as deprecated via 3
different mechanisms:
- Marking them with a deprecation warning `MarkParamaterDeprecatedWarning`. This
  will generate just a warning if the parameter is supplied..
- Marking them with a deprecation error `MarkParamaterDeprecatedError`. This
  will generate an error if the parameter is supplied..
- Marking them with as renamed `MarkParamaterRenamed`. This will generate an
  error if the parameter is supplied if the parameter is supplied.

\subsection DevManParametersSecFAQ6 FAQ6 Constrained ranges
Parameters in an input parameters block can be constrained to the following
- `chi_data_types::AllowableRangeList`, created with a `std::vector` of a
  specific type. The parameter can then only be one of those types.
- `chi_data_types::AllowableRangeLowLimit`. Parameter must be >= than the
  value specified when the range is specified as closed or just > that the value
  specified when the range is specified as open.
- `chi_data_types::AllowableRangeHighLimit` same as the low limit but just with
  a smaller than comparison.
- `chi_data_types::AllowableRangeLowHighLimit` combination of both a high and
  a low limit range.

Example:
\code
chi::InputParameters StrangeSolver::GetInputParameters()
{
  chi::InputParameters params =
    chi_physics::Solver::GetInputParameters();
  params.AddRequiredParameter<std::string>(
  "sdm_type", "The spatial discretization type to be used");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "sdm_type", AllowableRangeList::New({"FV", "PWLC", "PWLD"}));
  // etc.
\endcode
*/