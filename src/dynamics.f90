! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
module dynamics
    use dynamics_frequency_response
    use dynamics_rotation
    use dynamics_structural
    use dynamics_kinematics
    use dynamics_vibrations
    use dynamics_helper
    use dynamics_stability
    use dynamics_controls
    use dynamics_system_id
    use dynamics_linkage
    use dynamics_graph
    use dynamics_joints
    use dynamics_parallel_linkage
    use dynamics_geometry
    use dynamics_quaternions
    use dynamics_rigid_bodies
    use dynamics_line_elements
    use dynamics_variational_integrators
    use dynamics_linkage_dynamics
end module