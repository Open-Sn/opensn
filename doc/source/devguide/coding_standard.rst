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

Static constants should use upper case.

.. code-block:: c++

   const int MY_CONSTANT = 10;

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
