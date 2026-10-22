# Security Policy
This policy covers vulnerabilities in the code, workflows and Dockerfiles in this repository.


## Supported Versions
Only the latest `main` is supported. Tagged releases will be supported in the same way once they exist.


## Reporting a Vulnerability
Report vulnerabilities privately through GitHub private vulnerability reporting: open the repository's Security tab and choose "Report a vulnerability". Please do not open a public issue, pull request or discussion for a suspected vulnerability.

A report should say which file, workflow or image is affected, how to reproduce the problem, and what an attacker could do with it. The maintainers will acknowledge the report, work on a fix in private, and credit the reporter in the advisory unless asked not to.


## Scope
The workflows run in the user's own cloud project, under the user's own credentials, on the user's own data; this repository stores no credentials and ships no data beyond a small open-access test VCF. In scope are problems such as a workflow or script that exposes credentials, a Dockerfile that pulls an unverified artifact, or a CI workflow that can be abused to run untrusted code with write access. Vulnerabilities in third-party tools and base images should be reported to their own projects.
