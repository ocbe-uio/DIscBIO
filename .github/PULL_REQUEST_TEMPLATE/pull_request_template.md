---
name: Pull request
about: Create a pull request
title: ''
assignees: ''
---

Thank you for your contribution to DIscBIO! We kindly ask you to fill the Pull Request template below to streamline the review process.

**Warning:** This project follows [GitFlow](https://nvie.com/posts/a-successful-git-branching-model/), meaning a PR to master should always come from the latest `release-` branch. Please do not merge a feature branch directly into `master`, as this can cause asynchronies with the work of other contributors. You can choose the base of your PR in the buttons above.

# Summarize the changes introduced

A clear and concise description of what the bug is.

# Sanity checks

Please mark with an "x" if you can assert the following:
- [ ] The built package passes `devtools::check()` without errors or warnings
- [ ] The built package passes `BiocCheck::BiocCheck()` without errors or warnings

## If you haven't marked all the options above, please paste the errors/warnings found

```
Paste output here
```

# Additional information

Anything else you would like to add?