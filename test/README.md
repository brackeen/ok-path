##  Test (Linux, Mac)
```
mkdir tmp && cd tmp && cmake .. && cmake --build . && ctest --verbose; cd .. && rm -Rf tmp
```

## Generate project

NOTE: The `build` folder is ignored via `.gitignore`.

```
mkdir -p build && cd build && cmake -G Xcode ..
```
