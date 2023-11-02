/**\page DevManOOPvsFunc OOP vs Functional Programming
 *

## Object Oriented Programming or Functional Programming

### <CENTER>Chi-Tech uses a healthy balance between OOP and FP</CENTER>

The fundamental language with which a human provides instructions to a
 computer has always been **assembly language**. This language,
 which was pioneered halfway through the 20th century, is specific to
 a computer architecture and became the lowest level of programming that
 rose above simply supplying bytes of binary code to the computer.

Assembly language comprises combinations of keywords and/or numbers to
 form instructions, with each instruction telling the processor what
 action to take, whether it is to move data from one memory location
 to another or to perform an arithmetic operation on some data. By
 stringing together numerous instructions the complexity of a human
 thought-process could be made repeatable by a computer, however, the
 process of generating assembly language programs quickly became tedious.
 It did not take long before mankind comprehended the need for "packaging"
 multiple instructions into **macros**, which could be reused and
 hence the modern programming language was born where macros are
 standardized and collected into higher level languages
 (i.e., FORTRAN, C, etc).

In the 1950s scientific programming revolved around **block programming**,
 a programming paradigm where
 thought-processes are arranged as blocks that creates-,
 transforms- or uses data. This era of programming also included
 the use of the famous "GOTO" statement which led to the concept
 of what we today call "spaghetti-code", a style of programming
 only seen in older programs. The biggest concern with this style
 is the amount of effort required by a person, other than the creator
 of the programmer, to understand the overall program structure and flow.
 The difficulties are mainly associated with the limitations of the
 human mind to maintain scope and in many cases programmers would find
 themselves having trouble comprehending programs that they wrote a
 few weeks or months ago. These problems were greatly reduced with
 the advent of **structured programming** where the flow of the
 program is controlled by means of control blocks like if/then/else
 statements within which blocks of code, called "functions", could
 be executed. This aligned well with the industrial era where procedures
 could be regarded as stand-alone parts and the overall flow of the program
 could be visualized by a process diagram. The general term used
 today is **functional programming** (FP) and is very effective
 for certain algorithms.


The general motivation for functional programming stemmed from a
 development point of view. In the true sense of talking to the
 computer, the FP paradigm still resulted in sets of instructions
 to the computer, however, it was found that reusing, changing or
 understanding code using block programming proved to be a tremendously
 difficult task, especially when "GOTO" statements constantly interrupted
 the linear flow of a program. Functional programming overcame this
 difficulty by defining stand-alone procedures that had predefined inputs
 and outputs along with a concept called **local scope**, which is
 a tremendous improvement that allowed much flexibility with variable
 names and allowed the programmer to remove the dependence on global data.
 This paradigm is still very much valid today and greatly improved the
 process of code development.


Another programming paradigm, called **object oriented programming**
 (OOP), was created in the 1960s to allow programmers to structure
 their code as objects that interact with each other. It ushered-in
 an era where the complexity of a computer programming could be limitless
 and the ease of comprehension of computer programs could be increased
 by relating code concepts to real world objects. This paradigm grew
 very popular and successful. It allowed for an alternative to the
 control blocks used in procedural programming and allowed programs to
 reduce or even eliminate cache-misses (a performance aspect). It's
 success however, was not long-lived, as the difficulties associated
 with this paradigm became apparent. One such difficulty was that a
 program's evolution became dependent on the underlying architecture/layout
 of objects. If the architecture was flawed from the start the program's
 development would eventually grind to a halt, where additions or
 modifications became near impossible.

OOP allows for the implementation of advanced concepts some of them
plagued with tremendously risk. Concepts such as inheritance
and encapsulation offered many rewards whilst simultaneously introducing
problems.

<b>Inheritance</b> allows one to define a parent class of objects,
with associated procedures called methods, that provide basic
 functionality that can be used or specialized by child classes.
 A typical example is a parent class "solver" specialized
 into "diffusion-solver" or "transport-solver". The concept is
 useful in many circumstances and is generally very successful if
 not for the risk of creating categorical hierarchies that ultimately
 become unmanageable. Inheritance ultimately allows programmers
 to "program" their code into a corner, from which there is no escape
 other than to refactor large portions of
 architecture (if not the whole architecture). In scientific programs,
 where we need to solve a multitude of different kinds of problems,
 inheritance can create such a large hierarchical depth that the
 program becomes very dependent on the programmer's knowledge of
 the underlying architecture and less on the actual knowledge of
 problem-solving techniques. Therefore, when the program grows beyond
 a certain level of complexity, it can be more efficient to write
 specialized code from scratch than to devote the time to learn and
 comprehend the hierarchy. Inheritance is not dispensable however,
 and has a very unique place in modern programming. The ultimate
 solution is to keep scientific program to a very shallow level of
 inheritance in order to minimize the hierarchical depth of the
 architecture.

<b>Encapsulation</b> allows certain data items of an object to
 be **private** and therefore protected from incorrect
 manipulation. It has provided great utility and will be a key
 feature for many applications but unfortunately the concept also
 fuels the over-engineering of data structures. This is because
 objects that are over-encapsulated, cannot easily be modified by
 procedures living in a FP paradigm, reducing the compatibility of
 FP with OOP. Architectures that require perfect protection through
 encapsulation often exhibit an overwhelming focus on the actual
 "encapsulation" and less focus on the actual working code. There is
 a healthy trade-off between creating accessor/mutator
 functions (i.e., getters/setters) and how likely it is that unprotected
 access will cause undetected problems. In most cases scientific
 programmers should be more reliant on unit-tests as compared to
 encapsulation.
\n
\n
The modern focus of scientific programs should be a healthy balance
 between OOP and FP and therefore Chi-Tech follows such a paradigm.


 * */