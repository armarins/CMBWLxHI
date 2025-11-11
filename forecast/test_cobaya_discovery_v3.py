import cobaya
from cobaya.component import get_component_class, ComponentNotFoundError
import sys
import os
import importlib

print(f"Cobaya version: {cobaya.__version__}")
print(f"Cobaya installation path: {cobaya.__file__}")
print(f"Python executable used by this script: {sys.executable}")

print("\nPYTHONPATH environment variable:")
print(f"  Value: '{os.environ.get('PYTHONPATH')}'")

print("\nsys.path (Python's module search paths as seen by this script):")
# Corrected the sys.path printing loop
for i, p_path in enumerate(sys.path):
    print(f"  [{i}] {p_path}")

print("\nAttempting to get 'mcmc' sampler class (Cobaya's internal) via get_component_class:")
try:
    mcmc_class = get_component_class("mcmc", "sampler")
    print(f"  SUCCESS: Got mcmc class: {mcmc_class}")
except Exception as e:
    print(f"  FAILED to get mcmc class: {e}")

print("\nDirectly importing Cobaya's emcee wrapper module:")
try:
    # This is where Cobaya's EmceeSampler class is typically defined
    import cobaya.samplers.mcmc.emcee as cobaya_emcee_module
    print(f"  SUCCESS: Imported cobaya.samplers.mcmc.emcee: {cobaya_emcee_module}")
    if hasattr(cobaya_emcee_module, 'EmceeSampler'):
        print(f"    Found EmceeSampler class: {cobaya_emcee_module.EmceeSampler}")
    else:
        print(f"    WARNING: EmceeSampler class NOT found in the cobaya.samplers.mcmc.emcee module!")
except ImportError as e:
    print(f"  FAILED to import cobaya.samplers.mcmc.emcee: {e} (ImportError)")
except Exception as e_gen:
    print(f"  FAILED to import cobaya.samplers.mcmc.emcee: {e_gen} (General Error)")

print("\nAttempting to get 'emcee' sampler class via Cobaya's get_component_class (AGAIN, after attempting direct import):")
try:
    emcee_class_after = get_component_class("emcee", "sampler")
    print(f"  SUCCESS: Got emcee class: {emcee_class_after}")
except Exception as e:
    print(f"  FAILED to get emcee class: {e}")

print("\nDirectly importing Cobaya's ultranest wrapper module:")
try:
    # This is where Cobaya's UltraNestSampler class is typically defined
    import cobaya.samplers.nested.ultranest as cobaya_ultranest_module
    print(f"  SUCCESS: Imported cobaya.samplers.nested.ultranest: {cobaya_ultranest_module}")
    if hasattr(cobaya_ultranest_module, 'UltraNestSampler'):
        print(f"    Found UltraNestSampler class: {cobaya_ultranest_module.UltraNestSampler}")
    else:
        print(f"    WARNING: UltraNestSampler class NOT found in the cobaya.samplers.nested.ultranest module!")
except ImportError as e:
    print(f"  FAILED to import cobaya.samplers.nested.ultranest: {e} (ImportError)")
except Exception as e_gen:
    print(f"  FAILED to import cobaya.samplers.nested.ultranest: {e_gen} (General Error)")

print("\nAttempting to get 'ultranest' sampler class via Cobaya's get_component_class (AGAIN, after attempting direct import):")
try:
    ultranest_class_after = get_component_class("ultranest", "sampler")
    print(f"  SUCCESS: Got ultranest class: {ultranest_class_after}")
except Exception as e:
    print(f"  FAILED to get ultranest class: {e}")

print("\nDirectly importing base packages (confirmation):")
try:
    emcee_spec = importlib.util.find_spec("emcee")
    if emcee_spec:
        emcee_base = importlib.import_module("emcee")
        print(f"  SUCCESS: Imported base emcee package: version {emcee_base.__version__}, path {emcee_base.__file__}")
    else:
        print("  FAILED to find spec for base emcee package.")
except Exception as e:
    print(f"  Error importing base emcee: {e}")
try:
    ultranest_spec = importlib.util.find_spec("ultranest")
    if ultranest_spec:
        ultranest_base = importlib.import_module("ultranest")
        print(f"  SUCCESS: Imported base ultranest package: version {ultranest_base.__version__}, path {ultranest_base.__file__}")
    else:
        print("  FAILED to find spec for base ultranest package.")
except Exception as e:
    print(f"  Error importing base ultranest: {e}")