import pytest


@pytest.fixture(scope="function", autouse=True)
def set_up():
    import qcdb

    qcdb.driver.pe.clean_options()

    try:
        import psi4
    except ImportError:
        pass
    else:
        psi4.core.clean()
        psi4.core.clean_timers()
        psi4.core.clean_options()
        psi4.set_output_file("pytest_output.dat", True)
