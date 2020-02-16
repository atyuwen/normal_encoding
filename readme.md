Normal Vector Encoder/Decoder

---

- Encoding float3 normal vectors to uint16 by using optimized spherical coords together with alias method, with a bounded error less than 0.562432 degree.

- Both encoding and decoding are very fast. Encoding do use a search process to map jagged array entry to a linearized index, but the searching is short and cache friendly. (5 steps at most, <2 steps in average). Decoding don't do any searching at all.

- Encoder need a 3172-bytes look up table, while decoder only needs 1362 bytes. (can be further reduced to 1/2 if exploit symmetries)

- The code is not carefully optimized.

---

- Normal Encoding ref: J. Smith, G. Petrova, S. Schaefe, Encoding Normal Vectors using Optimized Spherical Coordinates.
http://faculty.cs.tamu.edu/schaefer/research/normalCompression.pdf

- Alias Method ref: Knuth, The Art of Computer Programming, Vol. 2, Sect. 3.4.1, ex. 7.

