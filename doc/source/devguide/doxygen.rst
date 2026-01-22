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

Coverage Requirements
---------------------

Mandatory documentation:

- All classes, structs, and unions.
- All public methods except for trivial getters/setters, such as:
  - Reference or constant/volatile reference to a member variable.
  - Methods that get or set the value of a class member.
- All member variables (private, protected, public).
- All ``static`` members and methods.
- All functions except for translation-unit local functions. Exceptions include:

  - ``static`` functions (C++98 approach):

    .. code-block:: cpp

       static void helpFunction() { ... }
  - Anonymous namespace functions (C++11 approach):

    .. code-block:: cpp

       namespace {
       void helpFunction() { ... }
       }

- All ``extern`` and ``constexpr`` variables declared in headers.
- All C++20 concepts.

Optional documentation:

- Protected methods: document if they are non-trivial.
- Private methods: document only if logic is complex or non-obvious.
- Deprecated APIs: only leave a brief note or no documentation at all.
- Enums: give a brief description of the enum purpose and a description for each enum value. Skip if
  the name is too obvious.
- Operator overloading: document if its intended usage is not self-explanatory.

Do not document:

- Files.
- Entities not declared in header files.
- Type aliasing (``using`` or ``typedef``).
- Namespaces.
- Preprocessor macros (macro constants and macro functions).
- Default constructors.
- Copy/move constructors and copy/move assignments.
- Destructors.

Generic guidelines
------------------

- Use backslash-style commands (``\return``, ``\param``) instead of @-style commands.

- For single-line comments, use ``///``. For multi-line comments, use ``/** */``.

- The Doxygen commands ``\brief`` and ``\details`` should not be used. Doxygen is configured with
  ``JAVADOC_AUTOBRIEF=ON``, so the text up to the first full stop is automatically interpreted as
  the brief description. Any optional detailed description should begin on a new line following the
  brief description. For example:

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

Functions/methods
~~~~~~~~~~~~~~~~~

- Function or class method documentation must include a **one-line** summary describing what the
  function does, using the **verb in its base form** (without "s" or "es").

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

- Include ``\param``, ``\return`` (if not ``void`` and the returned value is not trivial from the
  function name), and ``\throw`` / ``\exception`` (if applicable).

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

- Describe what the type represents and state all parameter constraints (e.g., arithmetic, movable,
  comparable) in a **noun phrase**.

C++20 concept
~~~~~~~~~~~~~

- All concepts must be documented.
- Use a **noun phrase** describing the requirements and focus on what the concept guarantees.
- If the concept's behavior is complex, provide additional explanation in the detailed description.

Macro-dependent implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- For any class, function, or variable affected by macros, the documentation must state:

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
