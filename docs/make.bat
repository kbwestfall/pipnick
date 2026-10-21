@echo off
REM make.bat for Sphinx documentation

REM You can set these variables from the command line.
set SPHINXAPI=sphinx-apidoc
set SPHINXAPIOPT=members,private-members,undoc-members,show-inheritance
set SPHINXOPTS=-aE -w .\sphinx_warnings.out
set SPHINXBUILD=sphinx-build
set BUILDDIR=_build
set STATICDIR=_static
set DOCTREE=%BUILDDIR%\doctrees
set LOCALFILES=%BUILDDIR%\* api\*.rst sphinx_warnings.out

REM Internal variables.
set ALLSPHINXOPTS=-d %DOCTREE% %SPHINXOPTS% .

if "%1" == "clean" goto clean
if "%1" == "apirst" goto apirst
if "%1" == "htmlonly" goto htmlonly
if "%1" == "html" goto html

goto end

:clean
    echo Cleaning up...
    rmdir /S /Q %LOCALFILES%
    goto end

:apirst
    echo Generating API documentation...
    set SPHINX_APIDOC_OPTIONS=%SPHINXAPIOPT%
    %SPHINXAPI% --separate -o .\api ..\nickelpipeline ..\nickelpipeline\version.py
    goto end

:htmlonly
    echo Building HTML documentation...
    %SPHINXBUILD% -b html %ALLSPHINXOPTS% %BUILDDIR%\html
    goto end

:html
    call %0 apirst
    call %0 htmlonly
    goto end

:end
    echo Done.
