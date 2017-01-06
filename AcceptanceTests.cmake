# Set absolute tolerance to be used for testing
Set (abs_tol 5.0e-8)
Set (rel_tol 1.0e-13)

# Input:
#   - casename: with or without extension
#
Macro (add_acceptance_test casename)

  String (REGEX REPLACE "\\.[^.]*" "" basename "${casename}")

  Add_Test (NAME    ToF_accept_${casename}_all_steps
            COMMAND runAcceptanceTest
            "case=${OPM_DATA_ROOT}/flow_diagnostic_test/eclipse-simulation/${basename}"
            "ref-dir=${OPM_DATA_ROOT}/flow_diagnostic_test/fd-ref-data/${basename}"
            "atol=${abs_tol}" "rtol=${rel_tol}")

EndMacro (add_acceptance_test)

If (NOT TARGET test-suite)
  Add_Custom_Target (test-suite)
EndIf ()

# Acceptance tests

Add_Acceptance_Test (SIMPLE_2PH_W_FAULT_LGR)
