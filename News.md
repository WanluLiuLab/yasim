# News of YASIM

## 0.1.0

Release 0.1.0 `3d1774c0` was launched on 02/24/2022.

Updates:

- First usable version.

## 0.1.1

Updates:

- `commonutils` refactored.
  - `commonutils.importer._silent_tqdm` was split from `commonutils.importer.tqdm_importer`.
  - `fd` in `commonutils.io.ArchiveIO` changed to `_fd`.
  - Created package `commonutils.stdlib_helper`.
- `bioutils` changes:
  - FastQ reader added to `bioutils.io` and `bioutils.typing`.
  - `bioutils.feature` refactored to `bioutils.io.feature`, `bioutils.datastructure` and `nioutils.typing.feature`.
- Other formatting errors.
