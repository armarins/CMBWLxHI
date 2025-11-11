import cobaya
from cobaya.component import get_component_class, ComponentNotFoundError
import sys
import os
import importlib # For checking if base packages are importable

print(f"Cobaya version: {cobaya.__version__}")
print(f"Cobaya installation path: {cobaya.__file__}") # Path to cobaya/__init__.py
print(f"Python executable used by this script: {sys.executable}")

print("\nPYTHONPATH environment variable:")
print(f"  Value: '{os.environ.get('PYTHONPATH')}'")

print("\nsys.path (Python's module search paths as seen by this script):")
for i, p in enumerate(sys.path):
    print(f"  [{i}] {p}")

print("\nAttempting to get 'mcmc' sampler class (Cobaya's internal) via get_component_class:")
try:
    mcmc_class = get_component_class("mcmc", "sampler")
    print(f"  SUCCESS: Got mcmc class: {mcmc_class}")
    print(f"           Located at: {mcmc_class.__module__}.{mcmc_class.__name__}")
except ComponentNotFoundError as e:
    print(f"  FAILED to get mcmc class (ComponentNotFoundError): {e}")
except Exception as e_gen:
    print(f"  FAILED to get mcmc class (General Error): {e_gen}")

print("\nAttempting to get 'emcee' sampler class via Cobaya's get_component_class:")
try:
    emcee_class = get_component_class("emcee", "sampler")
    print(f"  SUCCESS: Got emcee class: {emcee_class}")
    print(f"           Located at: {emcee_class.__module__}.{emcee_class.__name__}")
except ComponentNotFoundError as e:
    print(f"  FAILED to get emcee class (ComponentNotFoundError): {e}")
except Exception as e_gen:
    print(f"  FAILED to get emcee class (General Error): {e_gen}")

print("\nAttempting to get 'ultranest' sampler class via Cobaya's get_component_class:")
try:
    ultranest_class = get_component_class("ultranest", "sampler")
    print(f"  SUCCESS: Got ultranest class: {ultranest_class}")
    print(f"           Located at: {ultranest_class.__module__}.{ultranest_class.__name__}")
except ComponentNotFoundError as e:
    print(f"  FAILED to get ultranest class (ComponentNotFoundError): {e}")
except Exception as e_gen:
    print(f"  FAILED to get ultranest class (General Error): {e_gen}")

print("\nDirectly importing base packages (confirmation):")
try:
    emcee_spec = importlib.util.find_spec("emcee")
    if emcee_spec:
        emcee_base = importlib.import_module("emcee")
        print(f"  SUCCESS: Imported base emcee package: version {emcee_base.__version__}, path {emcee_base.__file__}")
    else:
        print("  FAILED to find spec for base emcee package (is it installed?).")
except ImportError as e:
    print(f"  FAILED to import base emcee package: {e}")
except Exception as e_gen:
    print(f"  Error during base emcee package import/check: {e_gen}")

try:
    ultranest_spec = importlib.util.find_spec("ultranest")
    if ultranest_spec:
        ultranest_base = importlib.import_module("ultranest")
        print(f"  SUCCESS: Imported base ultranest package: version {ultranest_base.__version__}, path {ultranest_base.__file__}")
    else:
        print("  FAILED to find spec for base ultranest package (is it installed?).")
except ImportError as e:
    print(f"  FAILED to import base ultranest package: {e}")
except Exception as e_gen:
    print(f"  Error during base ultranest package import/check: {e_gen}")