# Coding Standards


This page describes the coding standard of OpenSn.

## File names

Directory and file names should use [snake](#naming-conventions) style.

```
some_directory/file_name1.ext
some_directory/file_name2.ext
some_other_directory/another_file
```

## C++ conventions


### Macros

Macro names should use [Pascal](#naming-conventions) case, macro parameters should
use [snake](#naming-conventions) style

```c++
#define MacroDefinition(macro_parameter)
```

### Namespaces

Namespace names should use [snake](#naming-conventions) style.

```C++
namespace ns_one {

...

}
```

### Enums

Enum names should use [Pascal](#naming-conventions) style, enum values upper case.

```c++
enum OurEnumType {
  ENUM_VALUE_1,
  ENUM_VALUE_2,
  ...
}
```

### Static constants

Static constants should use upper case.

```c++
const int MY_CONSTANT = 10;
```

### Classes and Structs

Class names should use [Pascal](#naming-conventions) style.
Member variables should use [snake](#naming-conventions) style.
Private and protected member variables should have trailing `_` (underscore).
Member functions should use [Pascal](#naming-conventions) style.
Member function parameters should use [snake](#naming-conventions) style.


```c++
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
```

### Getters and Setters

Getters should use the `Get` prefix and setters should use the `Set` prefix.

```c++
class MyClass
{
public:
  Type GetMember() { return member_; }
  void SetMember(Type type) { member_ = type; }

private:
  Type member_;
};
```

### Numbers

- Decimal numbers are written with both whole and fraction part.
  Examples: `1234.12`, `1.0`, `0.0`, `1.0e-15`.

### Boolean operators

Boolean operators `or`, `and` and `not` should be used instead of `||`, `&&` and `!`.

### Pointers

Shared pointers (`std::shared`) are preferred over raw pointers.
Exception to this rule is when the code interacts with a 3rd party library like PETSc where
shared pointers simply don't exist.

### Conditionals

A space should be used after the keyword in a conditional statement. There
is no space inside parentheses. The statement block should not be enclosed
in braces if the condition fits on a single line.

```c++
if (a == b)
  a = 0;
```

If the condition spans multiple lines, the statement block may be enclosed in
braces for clarity.

```c++
if (std::none_of(my_container.begin(),
                 my_container.end(),
                 [](int val) { return val == 0; }))
{
  return;
}
```

### Comments

In-code comments should use `//`.

```c++
// in-code comment
call();
```

For [doxygen](https://www.doxygen.nl/)-style comments, use `///` for single-line comments
and `/** */` for multi-line comments with the following formatting:

```c++
/**
 * Function description
 *
 * \param i Description of the first parameter
 * \return Description of the return value
 */
double SomeMemberFunction(int i)

/// Single-line comment
int my_member_variable;
```

### Include directives

Preprocessor `#include` directives should be ordered as follows:

1. Header for this compilation unit (.h file that corresponds to .cc/.cpp file)
1. Other OpenSn headers
1. Non-standard, non-system libraries
1. C++ headers
1. C headers

There should **not** be empty lines separating the groups.

Example:
```c++
#include "modules/cfem_diffusion/cfem_diffusion_solver.h"
#include "framework/data_types/varying.h"
#include "petsc.h"
#include <string>
#include <map>
```

## Command-line parameters

Command line parameters used by the OpenSn binary or any OpenSn script should use
[kebab](#naming-conventions) style.

## References

### Naming Conventions
- Snake style: `this_is_snake_style`
- Kebab style: `this-is-kebab-style`
- Pascal style: `ThisIsPascalStyle`
