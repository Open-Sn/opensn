Coding Standards
================

This page describes the coding standard of OpenSn.

File names
----------

Directory and file names should use `snake` style (see :ref:`naming-conventions`).

.. code-block:: text

   some_directory/file_name1.ext
   some_directory/file_name2.ext
   some_other_directory/another_file

C++ conventions
---------------

Macros
~~~~~~

Macro names should use `Pascal` style, macro parameters should use `snake` style (see :ref:`naming-conventions`).

.. code-block:: c++

   #define MacroDefinition(macro_parameter)

Namespaces
~~~~~~~~~~

The topmost namespace is ``opensn``. We allow exactly two levels of namespaces:

.. code-block:: c++

   namespace opensn {

   namespace solver_impl {
   ...
   } // solver_impl

   } // opensn

Do **not** introduce additional namespace levels. If you need further subdivision, refactor by:

- Splitting into separate libraries/modules.
- Moving related code into classes or subdirectories.

Namespace names should use `snake` style:

.. code-block:: c++

   namespace ns_one {
   ...
   } // ns_one

Use anonymous namespace for variables with internal linkage:

.. code-block:: c++

   namespace opensn {
   namespace {
   ...
   // internal variables, functions, etc.
   ...
   }
   } // opensn

Enums
~~~~~

Enum names should use `Pascal` style, enum values upper case.

.. code-block:: c++

   enum OurEnumType {
     ENUM_VALUE_1,
     ENUM_VALUE_2,
     ...
   }

Static constants
~~~~~~~~~~~~~~~~

Static constants in global scope should use upper case.

.. code-block:: c++

   const int MY_CONSTANT = 10;
   constexpr double MY_DOUBLE = 12.765;


Static constants in local scope should use `snake` case.

   .. code-block:: c++

      {
        const double local_tolerance = 1e-8;
        constexpr double another_local_tolerance = 1e-12;
        ...
      }

Classes and Structs
~~~~~~~~~~~~~~~~~~~

Class names should use `Pascal` style.
Member variables should use `snake` style.
Private and protected member variables should have trailing underscores (`_`).
Member functions should use `Pascal` style.
Member function parameters should use `snake` style.

.. code-block:: c++

   class ThisIsAClassName {
   public:
      int public_member_var;
   protected:
      int my_member_variable_;
      double another_member_variable_;

      void MyCoolMemberFunction();
      void MemberFunctionWithAnArgument(int argument_name);
   };

   struct ThisIsAStructName {
      int public_member_var;

      void MyCoolMemberFunction();
      void MemberFunctionWithAnArgument(int argument_name);

   protected:
      int my_member_variable_;
      double another_member_variable_;
   };


The order of variables and functions inside a `class`/`struct`  should be as shown below:

.. code-block:: c++

   class ClassName {
   public:
      # member functions
      # member variables
   protected:
      # member functions
      # member variables
   private:
      # member functions
      # member variables

   public:
      # static function
      # static variables
   protected:
      # static function
      # static variables
   private:
      # static function
      # static variables
   };

Note: The is the *preffered* order. It is not always possible to achieve this in cases where
structs and enums must be declared before used. Those cases are allowed exceptions for
deviating from this ordering.

Getters and Setters
~~~~~~~~~~~~~~~~~~~

Getters should use the `Get` prefix and setters should use the `Set` prefix.

.. code-block:: c++

   class MyClass
   {
   public:
     Type GetMember() { return member_; }
     void SetMember(Type type) { member_ = type; }

   private:
     Type member_;
   };

Numbers
~~~~~~~

- Decimal numbers are written with both whole and fraction part.
  Examples: ``1234.12``, ``1.0``, ``0.0``, ``1.0e-15``.

Boolean operators
~~~~~~~~~~~~~~~~~

Boolean operators ``or``, ``and`` and ``not`` should be used instead of ``||``, ``&&`` and ``!``.

Pointers
~~~~~~~~

Shared pointers (``std::shared_ptr``) are preferred over raw pointers.
Exception to this rule is when the code interacts with a 3rd party library like PETSc where shared pointers simply don't exist.

Conditionals
~~~~~~~~~~~~

A space should be used after the keyword in a conditional statement.
There is no space inside parentheses. The statement block should not be enclosed in braces if the condition fits on a single line.

.. code-block:: c++

   if (a == b)
     a = 0;

If the condition spans multiple lines, the statement block may be enclosed in braces for clarity.

.. code-block:: c++

   if (std::none_of(my_container.begin(),
                    my_container.end(),
                    [](int val) { return val == 0; }))
   {
     return;
   }

Comments
~~~~~~~~

In-code comments should use ``//``.

.. code-block:: c++

   // in-code comment
   call();

For `doxygen <https://www.doxygen.nl/>`_-style comments, refer to :ref:`doxygen-guidelines` section.

Include directives
~~~~~~~~~~~~~~~~~~

Preprocessor ``#include`` directives should be ordered as follows:

1. Header for this compilation unit (.h file that corresponds to .cc/.cpp file)
2. Other OpenSn headers
3. Non-standard, non-system libraries
4. C++ headers
5. C headers

There should **not** be empty lines separating the groups.

Example:

.. code-block:: c++

   #include "modules/cfem_diffusion/cfem_diffusion_solver.h"
   #include "framework/data_types/varying.h"
   #include "petsc.h"
   #include <string>
   #include <map>

Lambdas
~~~~~~~

Named lambdas should use `Pascal` style.

.. code-block:: c++

   auto MyLambda = [](...) { ... };

Problem and solver validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These conventions apply to LBS problems and the solvers that drive them. Groupsets, sources, cross
sections, and other independent objects validate their own input.

Kinds of validation
^^^^^^^^^^^^^^^^^^^

A *problem rule* defines a valid combination of problem options and state. Put these rules in
``CheckConfigurationErrors``. An override calls the base implementation, then appends one message
for each additional violation. ``ValidateConfiguration`` runs all rules and throws one
``std::invalid_argument`` with the complete list.

A *solver requirement* defines which valid problems a solver can drive. Put discrete-ordinates
solver requirements in ``CheckRequirements``. For example, a k-eigenvalue solver requires
steady-state mode, no external source, and a fissionable material. Another solver may accept the
same problem, so these checks do not belong to the problem.

``Check...`` functions append errors without throwing. ``Validate...`` functions run checks and
throw. Keep these checks at their point of use:

* Report parse failures where input is read. Examples include a missing parameter or unknown
  boundary name.
* Validate a single argument where it is used. For example, ``SetTimeStep`` requires ``dt > 0``.
* Report failures found during work at that point. Examples include I/O errors and non-finite
  eigenvalues.

Constructors enforce the input schema while parsing. Checks that depend on the meaning or
combination of settings, the mesh, or the solver run after construction, when virtual dispatch is
safe.

Problem construction
^^^^^^^^^^^^^^^^^^^^

A problem is created through its static ``Create`` factory. The constructor parses and stores the
configuration. ``LBSProblem::Build`` then:

#. It validates the parsed configuration.
#. It builds common and derived runtime data.
#. It validates the built state.
#. It marks the problem built after both validation passes succeed.

A new problem class makes its constructor non-public and implements ``Create`` as follows:

.. code-block:: c++

   std::shared_ptr<MyProblem>
   MyProblem::Create(const ParameterBlock& params)
   {
     return Build(std::shared_ptr<MyProblem>(
       new MyProblem(MakeInputParameters<MyProblem>("lbs::MyProblem", params))));
   }

Put runtime setup in ``InitializeSpatialDiscretization`` and ``BuildRuntimeData``, not in the
constructor. A ``BuildRuntimeData`` override must call its base implementation.

Solver lifecycle
^^^^^^^^^^^^^^^^

``Solver::Initialize``, ``Execute``, and ``Advance`` are non-virtual wrappers. Each calls
``ValidateState`` before its derived-class hook. ``Initialize`` also validates after its hook
because initialization may change the problem. A failed initialization leaves the solver
uninitialized. ``Execute`` and ``Advance`` require a successful ``Initialize``.

``DiscreteOrdinatesSolver::ValidateState`` requires a built problem, runs its rules, then runs the
solver requirements. Other solver families implement ``ValidateState`` themselves. Validation
before every operation catches supported changes made after initialization.

Derived C++ solvers override ``InitializeSolver``, ``ExecuteSolver``, and ``AdvanceSolver`` rather
than the public lifecycle functions. They must also implement ``ValidateState``. Keep the public
wrappers non-virtual so derived solvers cannot bypass validation.

Runtime mutation
^^^^^^^^^^^^^^^^

A setter for a value read by problem rules first calls ``ValidateChange(member, value)``. This
helper validates the candidate, then restores the original member. A rule rejection therefore
leaves that member unchanged.

The setter then commits the value and calls ``RebuildRuntimeObjects``. That function refreshes
runtime objects shared by several setters.

.. code-block:: c++

   void DiscreteOrdinatesProblem::SetSaveAngularFlux(bool save)
   {
     ValidateChange(options_.save_angular_flux, save);
     options_.save_angular_flux = save;
     RebuildRuntimeObjects();
   }

A setter contains only value-specific work, such as remapping precursors after a cross-section
change. Source definitions need no runtime rebuild because source assembly reads them directly.

``ValidateChange`` is not a full transaction. It does not restore caller-owned objects or roll back
a later rebuild failure. Do not promise a stronger exception guarantee unless the setter provides
one.

Solvers sharing a problem
^^^^^^^^^^^^^^^^^^^^^^^^^

Several solvers may drive one problem in sequence. Flux, precursors, and time state remain in the
problem, so the next solver can continue from them.

A solver must restore temporary settings on problem-owned contexts. Use scoped guards such as
``WGSContext::OverrideSourceScopes``. Each solver also initializes its scratch data instead of
relying on a previous solver. Problem setters are persistent and affect every solver that shares
the problem.

MPI collectivity
^^^^^^^^^^^^^^^^

Problem construction, validation, supported setters, and discrete-ordinates solver calls are
collective. All ranks call them in the same order with equivalent input. A rule based on local data
must reduce it before deciding. Global cell counts include owned cells only, not ghost copies.


Command-line parameters
-----------------------

Command line parameters used by the OpenSn binary or any OpenSn script should use `kebab` style (see :ref:`naming-conventions`).

References
----------

.. _naming-conventions:

Naming Conventions
~~~~~~~~~~~~~~~~~~

- Snake style: ``this_is_snake_style``
- Kebab style: ``this-is-kebab-style``
- Pascal style: ``ThisIsPascalStyle``
