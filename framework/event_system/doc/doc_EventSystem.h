/**\defgroup doc_EventSystem I Event System

# ChiTech's event system
We have some elements in ChiTech that follow the `Publisher/Subscriber` design
pattern. The basic functionality is captured with the base classes
`chi::EventPublisher` and `chi::EventSubscriber`, where multiple subscribers are
assigned to a publisher.

\image html "EventSystem0.drawio.png" width=700px

## Events
Events are implemented as, a conceptually simple, data class `chi::Event` and
has the following members:
- `chi::Event::Name()`
- `chi::Event::Code()`
- `chi::Event::Parameters()`

The latter member is a `chi::ParameterBlock` supporting anything from simple
scalar quantities all the way through to complex hierarchies of data.

## System events
Systems events are posted via the `chi::SystemWideEventPublisher` singleton
and serves the purpose of being the central publisher where all events end up,
although, some publishers can attempt to handle events before this handler.
Events originating from this publisher are:
- `"ProgramStart"`
- `"ProgramExecuted"`

## Physics events
Physics events are handled via the `chi_physics::PhysicsEventPublisher`
singleton and has the following basic events (although this list will be
extended in future):
- `"SolverPreInitialize"`
- `"SolverInitialized"`
- `"SolverPreExecution"`
- `"SolverExecuted"`
- `"SolverPreStep"`
- `"SolverStep"`
- `"SolverPreAdvance"`
- `"SolverAdvanced"`

## `PhysicsEventPublisher` interaction with `SystemWideEventPublisher`
Order of operations and dependencies related to events can be resolved by using
different publishers. This is because the event publishers always form a
hierarchy with the `SystemWideEventPublisher` as the base. Any event published
by a leaf publisher will ultimately get forwarded to the
`SystemWideEventPublisher`. An example flow of events is shown below:

\image html "EventSystem1.drawio.png" width=700px


* */