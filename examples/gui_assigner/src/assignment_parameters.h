/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _ASSIGNMENT_PARAMETERS_H
#define _ASSIGNMENT_PARAMETERS_H

class AssignmentParameters {
public:
    double maxCloseEnoughRmsd;

    static auto getDefault() -> AssignmentParameters
    {
        return {
            .maxCloseEnoughRmsd = 0.5,
        };
    };
};


#endif // _ASSIGNMENT_PARAMETERS_H

