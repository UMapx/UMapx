# Source layout and documentation

- Keep new implementations and helper methods in the primary file of their class.
  Do not introduce `Class.Part.cs` files. The existing wavelet and Weyl–Heisenberg
  partial files are retained by the user's explicit preference.
- Document new numerical helpers with meaningful English XML comments: explain
  their purpose, parameters, result, and relevant domains or branch conventions.
  Explain non-obvious numerical steps near the code.
- Keep source code, comments, and project documentation in English; do not add
  Cyrillic characters.
