.. _dependency-maintenance:

Version and Dependency Maintenance
==================================

OpenSn keeps release metadata in a small set of authoritative files. This page
describes how those files are used and how to change them without leaving the
CMake build, Python package, or documentation out of sync.

Updating the OpenSn version
---------------------------

``VERSION.txt`` is the authoritative OpenSn version and must contain three
numeric components, ``MAJOR.MINOR.PATCH``, such as ``1.2.3``. CMake, the Python
package, and Sphinx read it directly, and CMake and ``setup.py`` reject any
other format. The build reconfigures automatically when the file changes.

The shared-library ``SOVERSION`` is the major version. Increment ``MAJOR`` for
changes that break the C++ or Python interface, ``MINOR`` for
backward-compatible features, and ``PATCH`` for backward-compatible fixes.

Making a release
----------------

A release is an annotated ``vMAJOR.MINOR.PATCH`` tag on a commit in ``main``,
plus a GitHub release for that tag. Publishing a release requires a maintainer
with permission to push tags to ``Open-Sn/opensn``. No workflow runs on tags or
releases, and the documentation site is deployed from ``main``, so publishing
the release does not start any builds.

The steps below use the remote names from :doc:`workflow`, where ``upstream``
is ``Open-Sn/opensn``.

#. **Open a release pull request.** On a branch from the latest ``main``:

   - Set ``VERSION.txt`` to the new version.
   - In ``distribution/spack/packages/opensn/package.py``, add
     ``version("X.Y.Z", tag="vX.Y.Z")`` above the previous release and point
     ``url`` at ``.../archive/refs/tags/vX.Y.Z.tar.gz``. Spack resolves the tag
     only when it fetches, so the entry can be added before the tag exists.
   - Update the ``opensn@`` examples in ``distribution/spack/README.md``.

   Merge the pull request after CI passes.

#. **Check the merged commit.** Fetch ``upstream`` and confirm that CI passed on
   the merge commit in ``main``. From a clean checkout of that commit, confirm
   that ``python setup.py --version`` and the CMake configure report the new
   version.

#. **Tag the release.** Create an annotated tag on the merge commit and push it:

   .. code-block:: shell

      git fetch upstream
      git tag -a vX.Y.Z -m "OpenSn vX.Y.Z" <merge-commit>
      git push upstream vX.Y.Z

   Do not move or delete a tag after it has been pushed. Anything that was
   fetched from it, including Spack, keeps the old commit. Fix a bad release
   with a new patch release.

#. **Publish the GitHub release.** With the GitHub CLI:

   .. code-block:: shell

      gh release create vX.Y.Z --repo Open-Sn/opensn --verify-tag \
        --title vX.Y.Z --generate-notes

   In the web interface, use **Releases** > **Draft a new release**, choose the
   existing tag, and select **Generate release notes**. Edit the generated
   notes into a short *What's Changed* list of user-visible changes and keep the
   *Full Changelog* comparison link. Do not upload files. GitHub attaches the
   source archives that the Spack ``url`` points to.
