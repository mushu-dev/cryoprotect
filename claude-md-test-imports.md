# Test CLAUDE.md Import Feature

This file tests the import syntax for CLAUDE.md files.

## Simple Imports

### API Module Import
@api/claude.md

### Database Module Import
@database/claude.md

## Contextual Notes

The imports above should be replaced with the content of their respective files when Claude Code processes this document.

## Nested Import Test

In theory, if one of the imported files also contains an import directive, it should be processed recursively.

## Path Formats

- Relative paths (as used above): `@api/claude.md`
- Absolute paths: `@/home/mushu/Projects/CryoProtect/chembl/claude.md`
- Dot notation: `@./pubchem/claude.md`

## Import Directive Variations

- Basic: `@path/to/file.md`
- With section: `@path/to/file.md#section`
- With line range: `@path/to/file.md[10:20]`