# Merges to dev
## Description
Provide a brief description of the PR's purpose here.

## Usage Changes
Document any changes to the parameters

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Pre-Review checklist (PR maker)
- [ ] Approval tests pass
- [ ] Acceptance tests pass
    - [ ] Running `setup.sh` and `run_test.sh` from acceptance_tests passes these tests:
      - [ ] All three modules run without errors 
      - [ ] ModCoor r^2 values converge to 1 with sufficient sampling (1000 samples)
      - [ ] ModCoor and CoorSpec agree (no subsampling). Based on plots.
    - [ ] New acceptance tests are documented 
- [ ] Linting isn't worse than it was (linting checks pass)
- [ ] Pytests pass
- [ ] New functions are documented
- [ ] New configuration parameters are documented (Usage changes above and added to Slack canvas)


## Review checklist (Reviewer)
- [ ] Approval tests pass
- [ ] Acceptance tests pass
    - [ ] Running `setup.sh` and `run_test.sh` from acceptance_tests passes these tests:
      - [ ] All three modules run without errors 
      - [ ] ModCoor r^2 values converge to 1 with sufficient sampling (1000 samples)
      - [ ] ModCoor and CoorSpec agree (no subsampling). Based on plots.
    - [ ] New acceptance tests are documented 
    - [ ] New acceptance tests are sufficient
    - [ ] My own test use-case(s) work(s)
- [ ] Linting scores haven't gotten worse (linting checks pass)
- [ ] Pytests pass: `bash run_pytests.sh`
- [ ] New functions are documented 
- [ ] New functions/variables are named appropriately
- [ ] No missed "low-hanging fruit" that would substantially aid readability.
- [ ] Any "high-hanging" or "rotten" fruit is added to the issues list.
- [ ] I understand what the changes are doing and how
- [ ] The PR itself is appropriately documented:
    - [ ] I understand the motivation for this PR
    - [ ] Usage changes are sufficiently documented

## Status
- [ ] Ready for review


# Merges to alpha
## Description
Provide a brief description of the PR's purpose here.

## Usage Changes
Document any changes to the parameters

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Pre-Review checklist (PR maker)
- [ ] Automated checks pass
- [ ] 0th order changes work. (Attach plots/outputs if appropriate)
- [ ] Pytests pass
- [ ] Changes are human-readable
- [ ] New configuration parameters are documented (Usage changes above)

## Review checklist (Reviewer)
- [ ] Automated checks pass
- [ ] I affirm that the 0th order changes work
- [ ] Pytests pass: `bash run_pytests.sh`
- [ ] I can read understand the changes
- [ ] I know what they're doing and how

## Status
- [ ] Ready for review
