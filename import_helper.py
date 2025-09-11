"""
Import helper to set up paths for pyanalib and other project modules
"""
import sys
import os

def setup_project_imports():
    """Add project root to Python path for imports using CAFPYANA_WD"""
    # Use CAFPYANA_WD environment variable if available
    cafpyana_wd = os.environ.get('CAFPYANA_WD')
    if cafpyana_wd:
        project_root = cafpyana_wd
    else:
        # Fallback to relative path from this file
        project_root = os.path.dirname(os.path.abspath(__file__))
    
    # Add to Python path if not already there
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
    return project_root

# Automatically set up imports when this module is imported
setup_project_imports()
