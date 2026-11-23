# Style Guide

This guide describes the conventions used in the `dynamics` Fortran library.
Follow the patterns in the surrounding code when they differ from this guide,
and prefer consistency within a module over sweeping reformatting. The rules
below reflect what the existing sources do; where the sources are inconsistent,
the guide states the preferred form for new code.

## Files and Naming

- Source files live under `src/`, grouped by subsystem (for example
  `src/rigid_bodies/`, `src/structural/`, `src/c/`).
- Name files and modules `dynamics_<topic>.f90` and use the same name for the
  module, for example `dynamics_flexible_assembly`.
- Test modules live in `test/` and use the `dynamics_<topic>_tests.f90` form.
- Begin new source files with the same copyright and SPDX license header used
  by existing files (see [dynamics_controls.f90](src/controls/dynamics_controls.f90)).

## Layout and Formatting

- Use free-form Fortran 2018 with the `.f90` extension.
- Indent with four spaces. Do not use tab characters; two existing files
  contain tabs and should not be used as a model.
- Keep lines at or below 132 characters.
- Write Fortran keywords in lowercase.
- Use a separator line of 80 characters (`! ---...`) between procedures
  after `contains`.
- Use `end subroutine` and `end function` without the procedure name, as in
  [dynamics_flexible_assembly.f90](src/flexible_bodies/dynamics_flexible_assembly.f90).

## Modules and Scoping

- Every module should contain `implicit none`.
- Make entities `private` by default and declare public names explicitly with
  `public ::`.
- Use `use <module>, only : <names>` rather than importing whole modules.
- Take kinds from `iso_fortran_env`, such as `real64` and `int64`. Do not use
  legacy forms such as `real(8)` or `real*8`.
- Do not use legacy constructs such as `goto`, `common`, `equivalence`, or the
  old operators `.eq.`, `.ne.`, `.lt.`, and so on. Use `==`, `/=`, `<`, and
  related symbols.

## Declarations

- Declare one entity per line for function and subroutine arguments:
  ```fortran
  real(real64), intent(in) :: x
  real(real64), intent(in) :: y
  ```
  Do not write `real(real64), intent(in) :: x, y`.
- Declarations can be combined for variables within a function or subroutine:
  `real(real64) :: x, y` is ok within the routine.
- For deferred-shape arrays, put bounds in the `dimension` attribute, not
  after the entity name. This applies to allocatable arrays and assumed-shape
  dummy arguments:
  ```fortran
  real(real64), intent(in), dimension(:) :: x
  real(real64), allocatable, dimension(:,:) :: y
  ```
  Fixed-size arrays may put bounds after the entity name:
  `real(real64) :: x(10)`.
- Give every dummy argument an explicit `intent`.
- Type-bound procedures take the object as `class(...), intent(in) :: this`
  (or `intent(inout)`, as needed).

## Documentation

Documentation comments use FORD-style `!!` markers. Plain `!` is for
implementation notes.

- Place a module description immediately after the `module` line.
- Place a procedure description as the first `!!` block inside the procedure,
  before the declarations.
- Place argument documentation as a `!!` block directly below each argument
  declaration. Describe the contract: units, allowed ranges, allocation state,
  and what the routine returns or rejects.
- Document behavior on invalid input in the procedure description, for example
  whether the routine returns a status or calls `error stop`.
- Only comment code where a reader would otherwise need clarification.
  Do not restate what the code obviously does.

## Errors and Invariants

- Use `error stop "<message>"` for violated invariants and invalid
  configurations. This is the convention used throughout the library.
- Do not use bare `stop`.
- Name the error condition clearly in the message, and document it in the
  procedure description.

## C Interface

- C-facing modules are in `src/c/` and are named `dynamics_c_<topic>.f90`.
- Exported procedures use the `c_` prefix and `bind(C, name = "c_...")`, for
  example `c_evaluate_transfer_function`.
- Keep the C-facing signatures in sync with the Fortran implementation and
  with any documentation in [c_docs/](c_docs/).

## Tests

- Place tests in `test/` and add them to the existing module for the
  subsystem where possible.
- Name test procedures for the behavior they check, for example
  `invalid_flexible_input`.
- Add tests only for behavior a change introduces or fixes. Cover the edge
  cases that would fail without the change.

## Validation Before Submitting

Run the build and tests with one of the following:

```sh
fpm build
fpm test
```

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Known Inconsistencies

These exist in the current code and are not blockers for new work:

- Deferred-shape arrays declared on the entity (`real(real64) :: x(:)`)
  instead of with `dimension(:)` appear in 26 declarations, in
  [dynamics_c_frequency.f90](src/c/dynamics_c_frequency.f90),
  [dynamics_controls.f90](src/controls/dynamics_controls.f90),
  [dynamics_beam_elements.f90](src/structural/dynamics_beam_elements.f90), and
  [dynamics_structural.f90](src/structural/dynamics_structural.f90).
  Fixed-size entity-style arrays are allowed and are not counted here.
- Argument declarations with several entities on one line
  (`real(c_double), intent(in) :: a, b`) appear in 99 lines across 15 files.
  Most are in `src/c/`, for example
  [dynamics_c_variational.f90](src/c/dynamics_c_variational.f90) (31 lines) and
  [dynamics_c_parallel_linkage.f90](src/c/dynamics_c_parallel_linkage.f90) (12 lines).
  Others include
  [dynamics_linkage_dynamics.f90](src/rigid_bodies/dynamics_linkage_dynamics.f90)
  (14 lines) and
  [dynamics_truss_elements.f90](src/structural/dynamics_truss_elements.f90)
  (4 lines).
  Combined declarations of local variables are allowed by this guide and are
  not included in this count.
- Two files contain tab characters:
  [dynamics_linkage_dynamics.f90](src/rigid_bodies/dynamics_linkage_dynamics.f90)
  and
  [dynamics_variational_integrators.f90](src/solvers/dynamics_variational_integrators.f90).
- One line exceeds 132 characters in
  [dynamics_frequency_response.f90](src/frequency/dynamics_frequency_response.f90).

Fix these only when touching the affected code.
