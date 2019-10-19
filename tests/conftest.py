import pytest 

@pytest.fixture(scope="function", autouse=True)
def set_up():
    import qcdb
    qcdb.driver.pe.clean_options()
