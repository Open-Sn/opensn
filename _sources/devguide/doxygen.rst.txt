.. _doxygen-guidelines:

Doxygen guidelines
==================

Purpose
-------

Provide developer-facing documentation for OpenSn using Doxygen.

Documentation must enable:

- New developers to understand code structure and intent.
- Existing developers to navigate APIs easily.
- Focus is on public and developer-relevant internal APIs, not just end-user manuals.

Making and reviewing documentation changes
------------------------------------------

Use this section before adding or correcting documentation. It is deliberately written as a common
workflow so that documentation is consistent regardless of who makes the change.

1. Identify the declaration and classify it using the coverage requirements below. The exclusions
   override the mandatory and optional rules. Otherwise, mandatory rules override optional rules.
   In particular, document a static class member or method regardless of its access level.
2. Establish the declaration's intent from its declaration, corresponding definition, call sites,
   tests, and related class documentation as needed. Do not make claims that are unsupported by the
   source. Do not invent units, ownership, thread-safety, exception behavior, valid ranges,
   performance guarantees, or algorithmic details. Request review when the required meaning cannot
   be established safely.
3. Attach one documentation block using the placement and content rules below. Document overloads
   independently unless their semantics are identical.
4. When correcting existing documentation, retain correct meaning and replace only inaccurate,
   incomplete, stale, or nonconforming text. Remove stale parameter, template-parameter, return,
   and exception entries, and do not attach duplicate documentation to a declaration.
5. Review the result: every ``\param`` must name a declared parameter, every ``\tparam`` must name a
   declared template parameter, and every documentation block must remain attached to its intended
   declaration. Keep the change limited to documentation; do not modify declarations,
   implementations, preprocessor conditions, include directives, or unrelated formatting.

Coverage Requirements
---------------------

Mandatory documentation:

- All class, struct, union, and enum definitions.
- All public methods except for trivial getters/setters, such as:
  - Reference or constant/volatile reference to a member variable.
  - Methods that get or set the value of a class member.
- All non-default public constructors.
- All member variables (private, protected, public).
- All ``static`` members and methods.
- All namespace-scope functions declared in a header, including inline functions. Do not add a second
  block to an out-of-class definition in the same header when its declaration is already documented.
  Do document an inline definition when it is the only declaration in the header. Exclude
  translation-unit-local functions, including:

  - ``static`` functions (C++98 approach):

    .. code-block:: cpp

       static void helpFunction() { ... }
  - Anonymous namespace functions (C++11 approach):

    .. code-block:: cpp

       namespace {
       void helpFunction() { ... }
       }

- All ``extern`` and ``constexpr`` variables declared in headers, including namespace-scope
  ``inline constexpr`` variables.
- All C++20 concepts.

Optional documentation:

- Protected methods: document if they are non-trivial.
- Private methods: document only if logic is complex or non-obvious.
- Protected and private non-default constructors: follow the corresponding protected or private method
  rule.
- Deprecated APIs: only leave a brief note or no documentation at all.
- Enum values: document each value whose behavior, value, or distinction from the other values is not
  clear from its identifier. When uncertain, document the value.
- Operator overloads and conversion operators: document them when their intended usage, operand
  meaning, result, or side effect is not self-explanatory.

Do not document:

- Files.
- Entities not declared in header files.
- Type aliasing (``using`` or ``typedef``).
- Namespaces.
- Preprocessor macros (macro constants and macro functions).
- Default constructors.
- Copy/move constructors and copy/move assignments.
- Destructors.
- Lambdas and local classes or local variables inside a function body.
- A friend declaration when the same free function has a documented namespace-scope declaration. In
  that case, document the namespace-scope declaration only.

Generic guidelines
------------------

- Use backslash-style commands (``\return``, ``\param``) instead of @-style commands.

- Place a documentation block immediately before the declaration it describes. Use ``///`` for a
  single-line block and ``/** */`` for a multi-line block. The only permitted trailing documentation
  form is ``///<`` on an enum value.

- For a templated, attributed, or constrained declaration, place the documentation before the
  ``template`` declaration or the first attribute, not between declaration specifiers and the name.
  For a declaration that is also its only definition in a header, place the documentation immediately
  before that definition. Do not document a source-file definition separately.

- Use a multi-line block when the documentation contains a detailed description or Doxygen command.
  Put the brief on the first text line. Put detailed text on subsequent lines, separate prose from a
  command block with a blank line, and put each command on its own line. For example:

  .. code-block:: cpp

     /**
      * Compute a cell-average source value.
      *
      * \param cell Cell for which to evaluate the source.
      * \param time Evaluation time in seconds.
      * \return Source value in particles per cubic centimeter per second.
      */
     double ComputeSource(const Cell& cell, double time) const;

- The Doxygen commands ``\brief`` and ``\details`` should not be used. End every brief description
  with a full stop. The documentation build must enable ``JAVADOC_AUTOBRIEF=YES`` so that the text up
  to the first full stop is interpreted as the brief description. Any optional detailed description
  should begin on a new line following the brief description. For example:

  .. code-block:: cpp

     /**
      * Brief description.
      * Long and more descriptive comments.
      * ...
      */

- Write in clear, concise English. Avoid implementation details unless necessary for correct usage.

- Potential pitfalls (such as manual resource deallocation, possible segmentation faults, or
  infinite loops) caused by misuse must be documented using the ``\note`` (what the user should
  notice) or ``\warning`` (what to avoid) command.

- Use the canonical forms ``\param``, ``\tparam``, ``\return``, and ``\throw`` for parameter,
  template-parameter, return-value, and exception documentation. Do not use ``@``-style commands,
  ``\retval``, ``\exception``, ``\brief``, or ``\details``. Other valid Doxygen commands and
  supported mathematical markup may be used when needed to communicate the API accurately.

Classes/structs
~~~~~~~~~~~~~~~

- Class or struct documentation must include a **one-line** summary of the class purpose. The brief
  description must be a **noun phrase**, not a verb or a complete sentence.

- Avoid redundant comments, such as “Class that represents…”.

- The detailed description is optional. Full sentences can be used. When provided, it should
  explain:

  - Design rationale
  - Usage guidelines
  - Valid conditions (if relevant)
  - Ownership (if relevant)
  - Thread-safety (if relevant)
  - Applied design patterns (if relevant)
  - Interaction with other classes (if relevant)
  - Performance considerations (if relevant)

- Document the class **only at its definition**, not at forward declarations.

- All class members must be documented **inside the class body**, not outside. Do not document
  methods defined outside the class in the header.

- For C++17 class template argument deduction (CTAD):

  - Only document the class definition.
  - Deduction guides should be documented either in the class description or in the appropriate
    constructor.
  - If additional clarification is needed for deduction guides, use ``//`` comments only.

  Example:

  .. code-block:: cpp

     /// Non-owning view (brief description).
     template <class T>
     struct View {
       /// Data pointer.
       T* data;
       /// Size.
       std::size_t size;

       constexpr View() : data(nullptr), size(0) {}

       /// Constructor from range iterators.
       template <class It>
       View(It first, It last)
       {
         ...
       }
     };

     // Deduce T from the iterator's value type
     template <class It>
     View(It, It) -> View<typename std::iterator_traits<It>::value_type>;

Unions
~~~~~~

- Document the union **only at its definition**, not at forward declarations.

- Provide a **one-line noun phrase** brief description of the union's purpose at the top of the
  union body.

- Document all fields inside the union body using **one-line noun phrase** brief descriptions.

- If the union contains complex or nested types, provide additional explanation in inline comments
  (``//``) only as needed.

- Do not document methods outside the union body.

Enums
~~~~~

- Document an enum only at its definition, not at a forward declaration.

- Use a one-line noun phrase for the enum's purpose. Place the documentation immediately before the
  ``enum`` declaration.

- Describe every enum value whose meaning is not clear from its name. Use ``///<`` after a value on
  the same line, or place a ``///`` comment immediately before the value. Keep each description as a
  one-line noun phrase.

  Example:

  .. code-block:: cpp

     /// Angular quadrature type identifier.
     enum class AngularQuadratureType
     {
       INVALID = 0,
       PRODUCT_QUADRATURE = 1,
       SLDFE_SQ = 2,          ///< Self-lumped discrete finite-element quadrature.
       LEBEDEV_QUADRATURE = 3 ///< Lebedev spherical quadrature.
     };

Functions/methods
~~~~~~~~~~~~~~~~~

- Function or class method documentation must include a **one-line** summary describing what the
  function does, using the **verb in its base form** (without "s" or "es").

- For a non-default constructor, use a one-line noun phrase that identifies the construction source
  or initialization performed. Do not use ``\return`` for constructors.

- Treat a friend function as a namespace-scope function. If its friend declaration is the only
  declaration in the header, document that declaration; otherwise document its namespace-scope
  declaration only.

- Treat ``operator()`` as a method. Treat other operators and conversion operators as optional as
  specified above. For a documented operator, state the operation and the meaning of its operands and
  result when those are not evident from the operator itself.

- Do not document a defaulted, deleted, copy, move, or destructor declaration. This exclusion also
  applies when the declaration is written explicitly in a header. A delegating constructor follows
  the rule for the constructor it declares, not the delegated-to constructor.

- Document functions/methods **only at their declaration** in the header file, not at their
  definition in the source file.

- The optional detailed description should explain:

  - Purpose
  - Algorithm
  - Input and output behavior (if applicable)
  - Assumptions (if applicable)
  - Side effects (if applicable)
  - Performance notes (if applicable)
  - Preconditions (if applicable)
  - Exceptions (if applicable)

- Use full sentences when possible.

- Use one ``\param <name>`` entry for every named parameter of a documented function or constructor.
  Describe its meaning and, when relevant, direction, ownership, units, valid range, and relationship
  to other parameters. Do not document unnamed parameters. State a default argument only when it
  changes observable behavior; do not merely repeat its literal value.

- Use ``\return`` for every documented non-``void`` function that is not an exempt trivial getter.
  Describe the meaning of the returned value and, when relevant, its units, ownership, lifetime, or
  sentinel/error meaning. Do not use ``\return`` for constructors, destructors, or ``void``
  functions.

- Use ``\throw <exception-type>`` only for exceptions that are part of the documented interface.
  State the condition that throws each exception. Do not infer exception guarantees from an
  implementation detail.

- A public override whose behavior is unchanged from its documented base declaration needs no
  separate block; it inherits the base documentation. Document the override when it changes
  semantics, preconditions, side effects, complexity, or exception behavior.

Variables/class members
~~~~~~~~~~~~~~~~~~~~~~~

- Extern variables or class members must include a **one-line noun phrase** summary describing the
  variable's nature or purpose.

- Documentation must **not** describe entities solely by naming a mathematical symbol. For example,
  the following is unacceptable:

  .. code-block:: cpp

     /// Alpha.
     int alpha;
     /// Value of x in the algorithm.
     float x;

  Exceptions are allowed if the variables are class members explicitly defined in a mathematical
  formulation or algorithm provided in the class documentation.

- Documentation should describe **what the variable represents or how it is used**, not just restate
  its symbolic name.

- Apply the same placement and noun-phrase rules to required namespace-scope ``extern`` and
  ``constexpr`` variables, including ``inline constexpr`` variables. For a constant, document the
  physical meaning, units, or contract rather than its literal value alone.

Template parameters
~~~~~~~~~~~~~~~~~~~

- Document template parameters using ``\tparam`` when:

  - The template parameter has semantic meaning beyond a generic STL type (``typename T``).
  - The parameter imposes constraints (``requires``), expectations (``if constexpr``), or specific
    behavior.
  - The intended use is not obvious from standard STL conventions.
  - The parameter affects algorithm behavior, ownership, or performance.

- Documentation may be omitted when:

  - The template parameter is a generic typename.
  - Usage strictly follows well-known STL conventions.
  - No additional assumptions or constraints apply.

- Describe what the template parameter represents and state applicable type, value, template, or pack
  constraints (e.g., arithmetic, movable, comparable) in a **noun phrase**.

- Use the declared parameter name in each ``\tparam`` entry, including non-type parameters, template
  template parameters, and parameter packs. Add the entry to the documentation block for the
  declaration that introduces the template parameter.

C++20 concept
~~~~~~~~~~~~~

- All concepts must be documented.
- Use a **noun phrase** describing the requirements and focus on what the concept guarantees.
- State material ``requires`` constraints and guarantees in the detailed description when they are not
  evident from the concept name.

Macro-dependent implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- For any class, struct, union, enum, function, or variable affected by macros, the documentation
  must state:

  - Macro dependency.
  - Behavioral difference.
  - Availability.

  Example:

  .. code-block:: cpp

     /**
      * Brief description...
      * - If ``__ABC__`` is defined, ...
      * - Otherwise, ...
      */
     void processData();

- If an API only exists under a macro, the mandatory wording is ``Only available when `MACRO_NAME`
  is defined``.

- If the macro is not defined for Doxygen parsing and the API must still be documented, add ``||
  defined(DOXYGEN_SHOULD_SKIP_THIS)``.

  Example:

  .. code-block:: cpp

     #if defined(__ABC__) || defined(DOXYGEN_SHOULD_SKIP_THIS)
     /**
      * Brief description...
      * Only available when ``__ABC__`` is defined.
      */
     int a;
     #elif !defined(__ABC__) || defined(DOXYGEN_SHOULD_SKIP_THIS)
     /**
      * Brief description...
      * Only available when ``__ABC__`` is not defined.
      */
     unsigned int x;
     #endif
