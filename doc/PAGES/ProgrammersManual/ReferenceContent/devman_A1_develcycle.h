/**\page DevManDevCycle Development Cycle
 *
\section devman_sec0 Development  Cycle

The development process of ChiTech is taylored for a small team of developers
and we will strive to keep to this philosophy as far as possible.
 We aim to provide stable and reliable releases to users whilst simultaneously
 keeping the development process fluent.

\subsection devman_sec0_1 Package update cycle - Primary: August
ChiTech releases comprise a packaged archive that represents a reviewed
snapshot of the code as determined by the principle developers. Update cycles
coincide approximately with U.S. University semesters:\n
 - No later than end of August and no earlier than end of July.
 - No later than end of February and no earlier than beginning of January.

\subsection devman_sec0_2 Development update cycle - Every 2 weeks
We as a team are commiting to reviewing Pull Requests (PRs) at a 2 week
interval. This means that a collaborator should not have to wait more than
two weeks for their pull request to be reviewed. This does not apply to
changes requested on previous pull requests which themselves will be
considered to be a new 2 week process.
\n
The basic development process is shown below:

\image html DevelopmentCycle.png " " width=600px

The development master branch will only be merged from the master branch of
 contributor forks. We will not merge PRs from feature branches. Contributors
 should merge feature branches to their respective master branches before
 issueing pull requests.
\n\n

<B>Example Pull-Request:</B>\n
Monday January 9 a contributor makes a complex PR. The last dedicated PR-review
 committee sitting was Wednesday January 4. The next dedicated sitting is then
 Wednesday January 18 (2 weeks later) and therefore if the committee does not
 find time before that date then the contributor will at maximum have to wait
 till January 18. If the review fails or requests additional changes then the
 timer is reset and the maximum wait time is again two weeks.

<B>Most likely course of action:</B>
Most likely the review committee members will start to individually respond
 to the PR before the big review sitting. The contributor will see GitHub code
 conversations for committee members which can be resolved on the go as
 shown below.

\image html ReviewConversation.png " " width=600px


\subsection devman_sec0_4 Regression tests
Before any contributer makes a pull request an <B>\f$\alpha\f$-regression</B>
 test-suite must be run. This is a regression test suite that can be run on a
 number of local machine processes (typically on the contributor's personal
 machine). On a weekly basis we will execute a high performance
 <B>\f$\beta\f$-regression</B> test-suite. This test-suite comprises high
 process-count tests necessary to ensure that the scalability of core features
 remain unaffected.  If the \f$\beta\f$-regression fails, all pull requests that
 were merged after the last successful \f$\beta\f$-regression will be reverted and
 the contributors notified.

\subsection devman_sec0_5 Code review criterion

 -# Changes have a clear purpose and are useful.
 -# The code passes the \f$\alpha\f$-regression tests.
 -# If appropriate, test cases are added to regression or unit test suites.
 We will not force unit tests without good cause.
 -# No memory leaks (checked with valgrind).
 -# Conforms to the ChiTech style guide.
 -# No degradation of performance or greatly increased memory usage. This is
 not a hard rule – in certain circumstances, a performance loss might be
 acceptable if there are compelling reasons.
 -# New features/input are documented. <B>Complex mathematical procedures
 must be accompanied by a whitepaper</B>.
 -# No unnecessary external software dependencies are introduced. An extensive
 justification needs to be provided for more dependencies.


\section devman_sec0_6 Documentation quality
Our policy on documentation is split into two parts: Code annotation and
 Technical documentation. It is our mission to create a sustainable modern
 code for which good quality documentation is essential.

\subsection devman_sec0_7 Code annotation
Code contributions need to be annotated appropriately. There is a strong
 connection to code-styling but over-all contributors must ensure their code
 is comprehensible to the principal developers and future contributors. We are
 not going to enforce specific annotation styles but we do enforce that the
 <B>annotation must be consistent within large portions of code</B>.
\n
\n
Function/methods need to have `doxygen` documentation at the place of
 definition. This policy will be <B>strongly</B> enforced for input
 language (i.e. lua function calls) as this goes into the input reference
 manual directly.

\subsection devman_sec0_8 Technical documentation
The fundamental mission on technical documentation is to have contributors
 provide whitepapers. These are latex-documents that fully explain
 mathematical derivations and implementations and may include pseudo-code
 that can be reference to by code.
 *
 *
 *
 *
 */