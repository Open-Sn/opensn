# Binding C++ classes to Lua

The OpenSn package provides a set of Lua bindings for the C++ classes.
The bindings are generated using the [LuaBridge3](https://github.com/kunitoki/LuaBridge3) library
which is embedded in the OpenSn repository.

This guide explains how to bind C++ classes into Lua.

All binding code is located in `lua/lib/opensn.cc` inside `opensnlua::Console::Bind`.

## Namespaces

Outside of the global namespace we have several of our own that group together the code functionality:

- `log` for logging events,
- `aquad` for angular quadrature,
- `logvol` for logical volumes,
- `mesh` for mesh-related functionality,
- `fieldfunc` for field functions,
- `xs` for cross sections,
- `post` for post-processing,
- `solver` for solvers,
- `lbs` for linear Boltzmann solvers.

The namespaces should not be changed unless directed by developers.
All new bindings will be created inside the existing namespaces.

Example code for a `aquad` namespace binding:

```cpp
luabridge::getGlobalNamespace(L)
  .beginNamespace("aquad")
  ...
  .endNamespace();
```

## Classes

Let's show how to bind a new class on an example.
Say we have a new class `MyLovelyClass` that inherits from a `LogicalVolume` class.
First, find the `logvol` namespace.
You will see that the first class there is binding the `LogicalVolume` class.
This is also the parent class of our new class.
We add a new section for binding our new class (typically at the end of the namespace block before the `.endNamespace()` call).

```cpp
.deriveClass<MyLovelyClass, LogicalVolume>("MyLovelyClass")
.endClass()
.beginClass<std::shared_ptr<MyLovelyClass>>("MyLovelyClassPtr")
.endClass()
```

This is the basic binding for a class and should be self-explanatory.

The second part of this block is for binding a shared pointer to the class and is essential for the correct
operation of the scripts.

## Member functions

To bind a member function of a class, we use the following syntax:

```cpp
.deriveClass<MyLovelyClass, LogicalVolume>("MyLovelyClass")
.addFunction("myFunction", &MyLovelyClass::myFunction)
.endClass()
```

To bind a static member function of a class, use the following syntax:

```cpp
.deriveClass<MyLovelyClass, LogicalVolume>("MyLovelyClass")
.addStaticFunction("myStaticFunction", &MyLovelyClass::myStaticFunction)
.endClass()
```

If our class was using `InputParameters` for construction, this would be the syntax we would use for binding of the `Create` method.

## Classes as Parameters

If our class would be used inside a Lua table (this is how `InputParameters` are represented in Lua),
then we need to teach the parsing code how to construct such a class.

This is done in `lua/lib/parse_table.cc` in the `SetBlockParam` function.

We find the section for contructing logical volumes and add the following code:

```cpp
else if (cls_name == "MyLovelyClass")
  block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
```

## Member variables

Public member variables can be bound into Lua as properties.
Let's see an example how it is done on the `Vector3` class.

```cpp
luabridge::getGlobalNamespace(L)
  .beginClass<Vector3>("Vector3")
  ...
  .addProperty("x", &Vector3::x)
  .addProperty("y", &Vector3::y)
  .addProperty("z", &Vector3::z)
  ...
  .endClass();
```

## Enums

Lua does not have an enum type, so we need to bind enum values as variables.

For example, if we had an enum:

```cpp
enum MyEnum {
  VALUE1,
  VALUE2,
  VALUE3
};
```

The binding code would look like this:

```cpp
luabridge::getGlobalNamespace(L)
  .addVariable("VALUE1", MyEnum::VALUE1)
  .addVariable("VALUE2", MyEnum::VALUE2)
  .addVariable("VALUE3", MyEnum::VALUE3);
```

## More details

For more details see the documentation of [LuaBridge3](https://kunitoki.github.io/LuaBridge3/Manual#2---accessing-c-from-lua)
