# First approach

Implement the algorithm described [here](https://www.sciencedirect.com/science/article/abs/pii/S0360835207000678?via%3Dihub)
- The algorithm will be implemented in C++ to optimize and parallelize it as
  much as possible
- Bindings will be implemented to create a UI


## TODO

- [X] Have a working version in a single `cpp` file
- [ ] Parallelize the code (try using `<execution>`) and std algorithms
- [ ] Restructure the code in cpp modules
- [ ] Introduce cross platform building with cmake
- [ ] Define python bindings
- [ ] Implemente a simple web app using streamlit

## Parallel execution

To be able to compile (on Ubuntu 22.04) with support for parallel algorithms and
c++23 features such as the `print` library, I needed to install the
`libtbb-dev`` package and compile with:
```bash
g++ -std=c++23 -ltbb main.cpp
```
I used the gcc 14.1.0 version.
