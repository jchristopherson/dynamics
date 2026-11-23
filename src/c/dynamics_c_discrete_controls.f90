module dynamics_c_discrete_controls
    use iso_c_binding
    use ieee_arithmetic, only : ieee_is_finite
    use dynamics_discrete_controls
    use dynamics_controls, only : state_space
    use dynamics_controls_tuning, only : integral_state_feedback_controller, controller_optimizer
    use dynamics_c_types, only : c_state_space_model, c_transfer_function, &
        c_alloc_state_space_model, c_alloc_transfer_function
    use dynamics_c_errors, only : c_api_error
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR, DYN_NULL_POINTER_ERROR, &
        DYN_MEMORY_ERROR, DYN_ARRAY_SIZE_ERROR
    implicit none
    private

    type, bind(C) :: c_discrete_optimizer_options
        real(c_double) :: tolerance
        integer(c_int) :: max_evaluations
        real(c_double) :: control_weight
    end type

contains

    ! Validate the common discrete-time argument.
    logical function valid_sample_time(ts) result(ok)
        real(c_double), intent(in) :: ts

        ok = .not.c_api_error(.not.ieee_is_finite(ts) .or. ts <= 0, &
            DYN_INVALID_INPUT_ERROR, "Sample time must be finite and positive.")
    end function

    ! Load a C matrix, accepting a null address only for an empty shape.
    logical function read_matrix(ptr, nr, nc, a) result(ok)
        type(c_ptr), intent(in) :: ptr
        integer, intent(in) :: nr
        integer, intent(in) :: nc
        real(c_double), allocatable, intent(out), dimension(:,:) :: a
        real(c_double), pointer, dimension(:,:) :: values

        ok = .true.
        if (nr == 0 .or. nc == 0) then
            allocate(a(nr, nc))
            return
        end if
        ok = .false.
        if (c_api_error(.not.c_associated(ptr), DYN_NULL_POINTER_ERROR, "Missing matrix data.")) return
        call c_f_pointer(ptr, values, [nr, nc])
        if (c_api_error(.not.all(ieee_is_finite(values)), DYN_INVALID_INPUT_ERROR, "Nonfinite matrix data.")) return
        a = values
        ok = .true.
    end function

    ! Load a realization without dereferencing empty matrices.
    logical function read_model(model, ts, sys) result(ok)
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), intent(in) :: ts
        type(discrete_state_space), intent(out) :: sys
        integer :: n
        integer :: m
        integer :: p

        ok = .false.
        if (c_api_error(.not.present(model), DYN_NULL_POINTER_ERROR, "Missing state-space model.")) return
        n = model%dimension
        m = model%n_inputs
        p = model%n_outputs
        if (c_api_error(n < 0 .or. m < 1 .or. p < 1, DYN_INVALID_INPUT_ERROR, &
            "Invalid state-space dimensions.")) return
        if (.not.read_matrix(model%A, n, n, sys%A)) return
        if (.not.read_matrix(model%B, n, m, sys%B)) return
        if (.not.read_matrix(model%C, p, n, sys%C)) return
        if (.not.read_matrix(model%D, p, m, sys%D)) return
        sys%sample_time = ts
        ok = .true.
    end function

    ! Allocate and fill the C result only after a realization is available.
    integer(c_int) function write_model(a, b, c, d, model) result(status)
        real(c_double), intent(in), dimension(:,:) :: a
        real(c_double), intent(in), dimension(:,:) :: b
        real(c_double), intent(in), dimension(:,:) :: c
        real(c_double), intent(in), dimension(:,:) :: d
        type(c_state_space_model), intent(inout) :: model
        real(c_double), pointer, dimension(:,:) :: values
        integer :: n
        integer :: m
        integer :: p

        status = -1
        if (c_api_error(.not.(all(ieee_is_finite(a)) .and. all(ieee_is_finite(b)) .and. &
            all(ieee_is_finite(c)) .and. all(ieee_is_finite(d))), DYN_INVALID_INPUT_ERROR, &
            "State-space transform produced nonfinite values.")) return
        n = size(a, 1)
        m = size(b, 2)
        p = size(c, 1)
        if (c_api_error(c_alloc_state_space_model(n, m, p, model) /= 0, DYN_MEMORY_ERROR, &
            "Could not allocate state-space result.")) return
        if (n > 0) then
            call c_f_pointer(model%A, values, [n, n])
            values = a
            call c_f_pointer(model%B, values, [n, m])
            values = b
            call c_f_pointer(model%C, values, [p, n])
            values = c
        end if
        call c_f_pointer(model%D, values, [p, m])
        values = d
        status = 0
    end function

    ! Load polynomial coefficients in ascending powers of z**(-1).
    logical function read_tf(tf, ts, sys) result(ok)
        type(c_transfer_function), optional, intent(in) :: tf
        real(c_double), intent(in) :: ts
        type(discrete_transfer_function), intent(out) :: sys
        real(c_double), pointer, dimension(:) :: b
        real(c_double), pointer, dimension(:) :: a

        ok = .false.
        if (.not.valid_sample_time(ts)) return
        if (c_api_error(.not.present(tf), DYN_NULL_POINTER_ERROR, "Missing transfer function.")) return
        if (c_api_error(min(tf%numerator%order, tf%denominator%order) < 0 .or. &
            max(tf%numerator%order, tf%denominator%order) == huge(0_c_int), &
            DYN_INVALID_INPUT_ERROR, "Invalid transfer-function orders.")) return
        if (c_api_error(.not.c_associated(tf%numerator%coefficients) .or. &
            .not.c_associated(tf%denominator%coefficients), DYN_NULL_POINTER_ERROR, &
            "Missing transfer-function coefficients.")) return
        call c_f_pointer(tf%numerator%coefficients, b, [tf%numerator%order + 1])
        call c_f_pointer(tf%denominator%coefficients, a, [tf%denominator%order + 1])
        if (c_api_error(.not.all(ieee_is_finite(a)) .or. .not.all(ieee_is_finite(b)) .or. &
            all(a == 0), DYN_INVALID_INPUT_ERROR, "Invalid transfer-function coefficients.")) return
        sys = discrete_transfer_function(b, a, ts)
        ok = .true.
    end function

    function c_c2d(model, ts, method, output) result(status) bind(C, name = "c_c2d")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        integer(c_int), value :: method
        type(c_state_space_model), optional, intent(inout) :: output
        integer(c_int) :: status
        type(discrete_state_space) :: src
        type(discrete_state_space) :: dst

        status = -1
        if (.not.valid_conversion(ts, method, output)) return
        if (.not.read_model(model, ts, src)) return
        if (model%dimension == 0) then
            status = write_model(src%A, src%B, src%C, src%D, output)
            return
        end if
        if (c_api_error(.not.all(ieee_is_finite(sum(abs(src%A * ts), dim = 1))) .or. &
            .not.all(ieee_is_finite(sum(abs(src%B * ts), dim = 1))), DYN_INVALID_INPUT_ERROR, &
            "Scaled state-space matrices are outside the numerical range.")) return
        if (method == DISCRETIZE_ZOH) then
            ! Keep the Padé scaling factor representable in the matrix exponential.
            if (c_api_error(max(maxval(sum(abs(src%A * ts), dim = 1)), &
                maxval(sum(abs(src%B * ts), dim = 1))) >= scale(1.0_c_double, maxexponent(ts) - 3), &
                DYN_INVALID_INPUT_ERROR, "Zero-order hold transform exceeds the scaling range.")) return
        end if
        dst = c2d(state_space(src%A, src%B, src%C, src%D), ts, method)
        status = write_model(dst%A, dst%B, dst%C, dst%D, output)
    end function

    function c_d2c(model, ts, method, output) result(status) bind(C, name = "c_d2c")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        integer(c_int), value :: method
        type(c_state_space_model), optional, intent(inout) :: output
        integer(c_int) :: status
        type(discrete_state_space) :: src
        type(state_space) :: dst
        logical :: ok

        status = -1
        if (.not.valid_conversion(ts, method, output)) return
        if (.not.read_model(model, ts, src)) return
        if (size(src%A, 1) == 0) then
            status = write_model(src%A, src%B, src%C, src%D, output)
            return
        end if
        call try_d2c(src, method, dst, ok)
        if (c_api_error(.not.ok, DYN_INVALID_INPUT_ERROR, &
            "Discrete-to-continuous transform is singular, outside its domain, or did not converge.")) return
        status = write_model(dst%A, dst%B, dst%C, dst%D, output)
    end function

    ! Validate common conversion arguments before calling the Fortran solver.
    logical function valid_conversion(ts, method, output) result(ok)
        real(c_double), intent(in) :: ts
        integer(c_int), intent(in) :: method
        type(c_state_space_model), optional, intent(in) :: output

        ok = .false.
        if (.not.valid_sample_time(ts)) return
        if (c_api_error(.not.present(output), DYN_NULL_POINTER_ERROR, "Missing conversion output.")) return
        if (c_api_error(method /= DISCRETIZE_ZOH .and. method /= DISCRETIZE_TUSTIN, DYN_INVALID_INPUT_ERROR, &
            "Unknown discretization method.")) return
        ok = .true.
    end function

    function c_discrete_state_space_step(model, ts, u, x, next, y) result(status) &
            bind(C, name = "c_discrete_state_space_step")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        real(c_double), optional, intent(in), dimension(*) :: u
        real(c_double), optional, intent(in), dimension(*) :: x
        real(c_double), optional, intent(out), dimension(*) :: next
        real(c_double), optional, intent(out), dimension(*) :: y
        integer(c_int) :: status
        type(discrete_state_space) :: sys
        integer :: n
        integer :: m
        integer :: p

        status = -1
        if (.not.valid_sample_time(ts)) return
        if (.not.read_model(model, ts, sys)) return
        n = size(sys%A, 1)
        m = size(sys%B, 2)
        p = size(sys%C, 1)
        if (c_api_error(.not.present(u) .or. .not.present(y) .or. &
            (n > 0 .and. (.not.present(x) .or. .not.present(next))), DYN_NULL_POINTER_ERROR, &
            "Missing input, state, or output buffer.")) return
        if (n == 0) then
            y(1:p) = matmul(sys%D, u(1:m))
        else
            y(1:p) = sys%evaluate_output(u(1:m), x(1:n))
            next(1:n) = sys%evaluate_next_state(u(1:m), x(1:n))
        end if
        status = 0
    end function

    function c_discrete_state_space_roots(model, ts, zeros, capacity, roots, count) result(status) &
            bind(C, name = "c_discrete_state_space_roots")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        logical(c_bool), value :: zeros
        integer(c_int), value :: capacity
        complex(c_double_complex), optional, intent(out), dimension(*) :: roots
        integer(c_int), optional, intent(out) :: count
        integer(c_int) :: status
        type(discrete_state_space) :: sys
        complex(c_double_complex), allocatable, dimension(:) :: values

        status = -1
        if (.not.valid_sample_time(ts)) return
        if (.not.read_model(model, ts, sys)) return
        if (zeros) then
            if (c_api_error(size(sys%B, 2) /= size(sys%C, 1), DYN_INVALID_INPUT_ERROR, &
                "State-space zeros require equal input and output counts.")) return
            values = sys%zeros()
        else
            values = sys%poles()
        end if
        status = write_roots(values, capacity, roots, count)
    end function

    function c_discrete_transfer_function_roots(tf, ts, zeros, capacity, roots, count) result(status) &
            bind(C, name = "c_discrete_transfer_function_roots")
        type(c_transfer_function), optional, intent(in) :: tf
        real(c_double), value :: ts
        logical(c_bool), value :: zeros
        integer(c_int), value :: capacity
        complex(c_double_complex), optional, intent(out), dimension(*) :: roots
        integer(c_int), optional, intent(out) :: count
        integer(c_int) :: status
        type(discrete_transfer_function) :: sys
        complex(c_double_complex), allocatable, dimension(:) :: values

        status = -1
        if (.not.read_tf(tf, ts, sys)) return
        if (zeros) then
            values = sys%zeros()
        else
            values = sys%poles()
        end if
        status = write_roots(values, capacity, roots, count)
    end function

    ! Return only the count on a query or if the supplied buffer is too short.
    integer(c_int) function write_roots(values, capacity, roots, count) result(status)
        complex(c_double_complex), intent(in), dimension(:) :: values
        integer(c_int), intent(in) :: capacity
        complex(c_double_complex), optional, intent(out), dimension(*) :: roots
        integer(c_int), optional, intent(out) :: count

        status = -1
        if (c_api_error(.not.present(count), DYN_NULL_POINTER_ERROR, "Missing root count output.")) return
        count = size(values)
        if (capacity == 0 .and. .not.present(roots)) then
            status = 0
            return
        end if
        if (c_api_error(.not.present(roots), DYN_NULL_POINTER_ERROR, "Missing root buffer.")) return
        if (c_api_error(capacity < count, DYN_ARRAY_SIZE_ERROR, "Root buffer is too short.")) return
        roots(1:count) = values
        status = 0
    end function

    function c_evaluate_discrete_transfer_function(tf, ts, n, omega, response) result(status) &
            bind(C, name = "c_evaluate_discrete_transfer_function")
        type(c_transfer_function), optional, intent(in) :: tf
        real(c_double), value :: ts
        integer(c_int), value :: n
        real(c_double), optional, intent(in), dimension(*) :: omega
        complex(c_double_complex), optional, intent(out), dimension(*) :: response
        integer(c_int) :: status
        type(discrete_transfer_function) :: sys

        status = -1
        if (.not.read_tf(tf, ts, sys)) return
        if (c_api_error(n < 0, DYN_ARRAY_SIZE_ERROR, "Negative frequency count.")) return
        if (n > 0) then
            if (c_api_error(.not.present(omega) .or. .not.present(response), DYN_NULL_POINTER_ERROR, &
                "Missing frequency or response array.")) return
            if (c_api_error(.not.all(ieee_is_finite(omega(1:n) * ts)), DYN_INVALID_INPUT_ERROR, &
                "Nonfinite scaled frequency.")) return
            response(1:n) = sys%evaluate(omega(1:n))
        end if
        status = 0
    end function

    function c_discrete_state_space_to_transfer_function(model, ts, output) result(status) &
            bind(C, name = "c_discrete_state_space_to_transfer_function")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        type(c_transfer_function), optional, intent(inout) :: output
        integer(c_int) :: status
        type(discrete_state_space) :: sys
        type(discrete_transfer_function) :: tf
        real(c_double), pointer, dimension(:) :: values
        real(c_double), allocatable, dimension(:) :: a
        real(c_double), allocatable, dimension(:) :: b
        integer :: nb
        integer :: na

        status = -1
        if (.not.valid_sample_time(ts)) return
        if (c_api_error(.not.present(output), DYN_NULL_POINTER_ERROR, "Missing transfer-function output.")) return
        if (.not.read_model(model, ts, sys)) return
        if (c_api_error(size(sys%D, 1) /= 1 .or. size(sys%D, 2) /= 1, DYN_INVALID_INPUT_ERROR, &
            "Transfer-function conversion requires a SISO model.")) return
        tf = sys%to_transfer_function()
        b = tf%Y%get_all()
        a = tf%X%get_all()
        if (c_api_error(.not.all(ieee_is_finite(a)) .or. .not.all(ieee_is_finite(b)), &
            DYN_INVALID_INPUT_ERROR, "Transfer-function conversion produced nonfinite coefficients.")) return
        nb = size(b) - 1
        na = size(a) - 1
        if (c_api_error(c_alloc_transfer_function(nb, na, output) /= 0, DYN_MEMORY_ERROR, &
            "Could not allocate transfer-function result.")) return
        call c_f_pointer(output%numerator%coefficients, values, [nb + 1])
        values = b
        call c_f_pointer(output%denominator%coefficients, values, [na + 1])
        values = a
        status = 0
    end function

    function c_discrete_transfer_function_to_state_space(tf, ts, observable, output) result(status) &
            bind(C, name = "c_discrete_transfer_function_to_state_space")
        type(c_transfer_function), optional, intent(in) :: tf
        real(c_double), value :: ts
        logical(c_bool), value :: observable
        type(c_state_space_model), optional, intent(inout) :: output
        integer(c_int) :: status
        type(discrete_transfer_function) :: src
        type(discrete_state_space) :: dst

        status = -1
        if (c_api_error(.not.present(output), DYN_NULL_POINTER_ERROR, "Missing realization output.")) return
        if (.not.read_tf(tf, ts, src)) return
        if (c_api_error(src%X%get(1) == 0, DYN_INVALID_INPUT_ERROR, &
            "Canonical conversion requires a nonzero denominator constant.")) return
        if (observable) then
            dst = src%to_ocf_state_space()
        else
            dst = src%to_ccf_state_space()
        end if
        status = write_model(dst%A, dst%B, dst%C, dst%D, output)
    end function

    ! Load the discrete plant and column-major state/integral gains.
    logical function read_controller(model, ts, k, ki, plant, ctrl) result(ok)
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), intent(in) :: ts
        type(c_ptr), intent(in) :: k
        type(c_ptr), intent(in) :: ki
        type(discrete_state_space), intent(out) :: plant
        type(integral_state_feedback_controller), intent(out) :: ctrl
        integer :: m

        ok = .false.
        if (.not.valid_sample_time(ts)) return
        if (.not.read_model(model, ts, plant)) return
        m = size(plant%B, 2)
        if (.not.read_matrix(k, m, size(plant%A, 1), ctrl%K)) return
        if (.not.read_matrix(ki, m, size(plant%C, 1), ctrl%Ki)) return
        ok = .true.
    end function

    function c_discrete_integral_feedback_apply(model, ts, k, ki, output) result(status) &
            bind(C, name = "c_discrete_integral_feedback_apply")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        type(c_ptr), value :: k
        type(c_ptr), value :: ki
        type(c_state_space_model), optional, intent(inout) :: output
        integer(c_int) :: status
        type(discrete_state_space) :: plant
        type(discrete_state_space) :: loop
        type(integral_state_feedback_controller) :: ctrl

        status = -1
        if (c_api_error(.not.present(output), DYN_NULL_POINTER_ERROR, "Missing closed-loop output.")) return
        if (.not.read_controller(model, ts, k, ki, plant, ctrl)) return
        loop = ctrl%apply_to_discrete_state_space_plant(plant)
        status = write_model(loop%A, loop%B, loop%C, loop%D, output)
    end function

    function c_discrete_integral_feedback_analyze(model, ts, k, ki, stable, ise) result(status) &
            bind(C, name = "c_discrete_integral_feedback_analyze")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        type(c_ptr), value :: k
        type(c_ptr), value :: ki
        logical(c_bool), optional, intent(out) :: stable
        real(c_double), optional, intent(out) :: ise
        integer(c_int) :: status
        type(discrete_state_space) :: plant
        type(integral_state_feedback_controller) :: ctrl

        status = -1
        if (c_api_error(.not.present(stable) .and. .not.present(ise), DYN_NULL_POINTER_ERROR, &
            "Request stability, integral squared error, or both.")) return
        if (.not.read_controller(model, ts, k, ki, plant, ctrl)) return
        if (present(stable)) stable = ctrl%is_stable(plant)
        if (present(ise)) ise = ctrl%integral_squared_error(plant)
        status = 0
    end function

    function c_discrete_integral_feedback_optimize(model, ts, k, ki, options, cost) result(status) &
            bind(C, name = "c_discrete_integral_feedback_optimize")
        type(c_state_space_model), optional, intent(in) :: model
        real(c_double), value :: ts
        type(c_ptr), value :: k
        type(c_ptr), value :: ki
        type(c_discrete_optimizer_options), optional, intent(in) :: options
        real(c_double), optional, intent(out) :: cost
        integer(c_int) :: status
        type(discrete_state_space) :: plant
        type(integral_state_feedback_controller) :: ctrl
        type(controller_optimizer) :: opt
        real(c_double), pointer, dimension(:,:) :: gains
        real(c_double) :: value

        status = -1
        if (.not.read_controller(model, ts, k, ki, plant, ctrl)) return
        if (present(options)) then
            if (c_api_error(.not.ieee_is_finite(options%tolerance) .or. options%tolerance <= 0 .or. &
                options%max_evaluations < 1 .or. .not.ieee_is_finite(options%control_weight) .or. &
                options%control_weight < 0, DYN_INVALID_INPUT_ERROR, "Invalid controller optimization settings.")) return
            opt%tolerance = options%tolerance
            opt%max_evaluations = options%max_evaluations
            opt%control_weight = options%control_weight
        end if
        if (c_api_error(.not.ctrl%is_stable(plant), DYN_INVALID_INPUT_ERROR, &
            "Initial gains must stabilize the discrete plant.")) return
        call opt%optimize(plant, ctrl, value)
        if (c_api_error(.not.ieee_is_finite(value) .or. .not.ctrl%is_stable(plant), DYN_INVALID_INPUT_ERROR, &
            "Optimizer did not return a finite, stabilizing result.")) return
        if (size(ctrl%K) > 0) then
            call c_f_pointer(k, gains, shape(ctrl%K))
            gains = ctrl%K
        end if
        call c_f_pointer(ki, gains, shape(ctrl%Ki))
        gains = ctrl%Ki
        if (present(cost)) cost = value
        status = 0
    end function
end module
