import pytest 

@pytest.fixture(scope="function", autouse=True)
def set_up():
    import qcdb
    #qcdb.driver.pe.clean_options()
    qcdb.driver.pe.clean_nu_options()
#    psi4.set_output_file("pytest_output.dat", True)

