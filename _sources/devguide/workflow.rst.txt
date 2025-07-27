Developer Workflow
==================

To participate in the OpenSn development, you will need to fork the upstream
repository on GitHub.  Visit the
`OpenSn repo on GitHub <https://github.com/Open-Sn/opensn>`_ and
`create your personal fork <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_.

Then, clone the fork onto the machine where you will be developing:

.. code-block:: shell

   $ git clone git@github.com:<username>/opensn.git

where ``<username>`` is your GitHub username.  
(This assumes you have
`set up your SSH keys for accessing GitHub <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_.)

Your fork will be referred to as ``origin`` (both in the *git* commands and
in the text below).  Also set up an ``upstream`` remote for pulling changes
from the upstream repo you forked from:

.. code-block:: shell

   $ git remote add upstream git@github.com:Open-Sn/opensn.git

OpenSn development uses the pull-request workflow that is common on GitHub.
This means **all development is done on feature branches**.  
Do **not** develop on the ``main`` branch—you will most likely create problems
for yourself down the road.

.. admonition:: Attention

   Do **NOT** merge ``main`` into your feature branch.

Formating
---------

* **C++** code should be formatted according to our *clang-format* choices
  found in the
  `\.clang\_format <https://github.com/Open-Sn/opensn/blob/main/.clang-format>`_
  file at the repository root.

* **Python** linting and style-guide checking is performed with ``flake8``;
  our configuration lives in
  `\.flake8 <https://github.com/Open-Sn/opensn/blob/main/.flake8>`_
  at the repository root.

Create a Branch
---------------

First, make sure you have no local changes:

.. code-block:: shell

   $ git status
   On branch main
   Your branch is up to date with 'origin/main'.

   nothing to commit, working tree clean

Branch off of ``main``:

.. code-block:: shell

   $ git checkout -b <branch name> main
   Switched to a new branch '<branch name>'

Tips for naming your branch:

* Use lower-case letters.
* Dashes (``-``) work better than underscores (``_``).
* If working on an issue, you can name your branch ``issue/<issue number>``.

**Useful commands**

.. list-table::
   :header-rows: 1
   :widths: 35 25

   * - Action
     - git command
   * - Get a list of branches
     - ``git branch -a``
   * - See your active branch
     - ``git branch``

Create a Commit
---------------

In *git* you must **stage** your changes first.  
Then you commit those changes locally, creating a *commit* (also called a
*patch*).  Every commit has a SHA-1 hash that uniquely identifies it.

Example commit header:

.. code-block:: shell

   commit <SHA1>
   Author: A. U. Thor <a.u.thor@somewhere.com>
   Date:   <Date>

**Useful commands**

.. list-table::
   :header-rows: 1
   :widths: 55 25

   * - Action
     - git command
   * - Add (stage) a new file
     - ``git add <file>``
   * - Add modifications on a file
     - ``git add <file>``
   * - Add all new files and local changes
     - ``git add -A``
   * - Add only modified files
     - ``git add -u``
   * - Move a file
     - ``git mv <source> <destination>``
   * - Remove a file
     - ``git rm <file>``
   * - See status (staged, unstaged, untracked)
     - ``git status``
   * - Commit locally
     - ``git commit``

Commit-message template:

.. code-block:: text

   Short description (up to 70 chars)

   Detailed description.  List as much useful and relevant
   information as possible.  This can span multiple lines
   and even paragraphs.

Tips for commits:

* Keep commits small—easier and faster to review.
* Stick to **one topic per commit**.
* Avoid doing multiple things in one commit.
* Prefer a series of small patches over one large patch.

More handy commands:

.. list-table::
   :header-rows: 1
   :widths: 45 25

   * - Action
     - git command
   * - See what files have changed
     - ``git status``
   * - List commits
     - ``git log``
   * - Show code changes in commits
     - ``git log -p``
   * - Show local changes not yet staged
     - ``git diff``

Sending a Pull Request
----------------------

When you think you are finished with your branch, send your changes for
review.

Push the branch to your fork on GitHub:

.. code-block:: shell

   $ git push origin <branch name>
   Enumerating objects: 38, done.
   Counting objects: 100% (38/38), done.
   ...

#. Navigate to ``https://github.com/<username>/opensn.git``.  
#. `Create a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_
   targeting the ``main`` branch.

Someone will review your branch and may leave comments; otherwise, it will be
merged. You’ll receive an email notification either way.

If your branch is merged, update your ``main`` branch:

.. code-block:: shell

   $ git checkout main
   $ git pull upstream main

Then delete your local branch:

.. code-block:: shell

   $ git branch -d <branch name>

Fixing Your Branch
------------------

The easiest way to fix your branch is with *fixup* commits (see
`fun with autosquash <https://technosorcery.net/blog/2010/02/fun-with-the-upcoming-1-7-release-of-git-rebase---interactive---autosquash/>`_):

#. Go patch-by-patch and fix code according to reviewer comments.  
#. Stage your changes (typically ``git add -u``).  
#. Commit with ``git commit --fixup=<commit>`` where ``<commit>`` is the
   SHA-1 of the patch you are fixing.  
#. Repeat until everything is fixed.  
#. Run ``git rebase --interactive --autosquash main``.  Save and quit the
   editor; *git* will do the rest.  
#. Verify with ``git log -p``.  
#. Push the updated branch: ``git push -f``

**Useful commands**

.. list-table::
   :header-rows: 1
   :widths: 50 25

   * - Action
     - git command
   * - Create a fixup patch
     - ``git commit --fixup=<sha1>``
   * - Add staged changes to the top-most commit
     - ``git commit --amend --no-edit``

Updating Your Branch
--------------------

If you need the latest ``main`` (and it has changed since you started your
branch), rebase your branch on top of it:

.. code-block:: shell

   $ git checkout main
   $ git pull upstream main
   $ git checkout <my-branch-name>
   $ git rebase main

*Note*: You will likely encounter conflicts that you must resolve.

Squashing Commits
-----------------

If asked to squash your commits:

.. code-block:: shell

   $ git rebase -i main

Your editor opens with a list of commits—change ``pick`` to ``f`` (fixup) for
all but the first commit, e.g.:

.. code-block:: text

   pick <sha1> Short description
   f <sha1> Some description
   f <sha1> Some other description

Save and quit; *git* squashes the commits into one.

You may edit the commit message with ``git commit --amend`` and then push:

.. code-block:: shell

   $ git push -f

.. admonition:: Attention

   This assumes you followed this guide and worked in your branch
   ``branch name`` and **not** in ``main``.
