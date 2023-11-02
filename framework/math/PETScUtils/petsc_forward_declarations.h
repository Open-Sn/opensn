#ifndef CHITECH_PETSC_FORWARD_DECLARATIONS_H
#define CHITECH_PETSC_FORWARD_DECLARATIONS_H

/**\file petsc_forward_declarations.h
 * This header file provides convenient forward declarations for common
 * PETSc types. This negates the need to include the bulky PETSc headers.*/

// NOLINTBEGIN(bugprone-reserved-identifier)
typedef struct _p_Vec* Vec;
typedef struct _p_Mat* Mat;
// NOLINTEND(bugprone-reserved-identifier)

#endif // CHITECH_PETSC_FORWARD_DECLARATIONS_H
