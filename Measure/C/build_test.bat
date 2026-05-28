@echo off
REM ========================================
REM Build script for Mirror Reflection 2D Simulation
REM ========================================

echo ========================================
echo Building Mirror Reflection 2D Simulation
echo ========================================

REM Check if MSVC compiler is available
where cl >nul 2>&1
if %errorlevel% equ 0 (
    echo Using MSVC compiler...
    echo
    
    REM Compile with optimization
    cl /O2 /W3 test_mirror_2d.c mirror_reflection_simulation_2d.c /Fe:test_mirror_2d.exe /link
    
    if %errorlevel% equ 0 (
        echo
        echo ========================================
        echo Build successful!
        echo ========================================
        echo
        echo Running test...
        echo
        test_mirror_2d.exe
    ) else (
        echo
        echo Build failed!
        exit /b 1
    )
) else (
    echo MSVC compiler not found.
    echo Please install Visual Studio or use GCC.
    echo
    
    REM Try GCC
    where gcc >nul 2>&1
    if %errorlevel% equ 0 (
        echo Using GCC compiler...
        echo
        
        gcc -O2 -Wall -o test_mirror_2d.exe test_mirror_2d.c mirror_reflection_simulation_2d.c -lm
        
        if %errorlevel% equ 0 (
            echo
            echo ========================================
            echo Build successful!
            echo ========================================
            echo
            echo Running test...
            echo
            test_mirror_2d.exe
        ) else (
            echo
            echo Build failed!
            exit /b 1
        )
    ) else (
        echo GCC compiler not found.
        echo Please install a C compiler (MSVC or GCC).
        exit /b 1
    )
)

echo
echo ========================================
echo Done!
echo ========================================
