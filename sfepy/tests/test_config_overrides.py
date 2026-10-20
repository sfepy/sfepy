from copy import deepcopy

from sfepy.base.base import update_dict_recursively
from sfepy.base.conf import ProblemConf


def test_override_field_tuple():
    data = {'fields': {'temperature': ('real', 'scalar', 'Omega', 1)}}
    override = {'fields': {'temperature': ('real', 'scalar', 'Omega', 3)}}
    conf = ProblemConf.from_dict(data, None, required=[], override=override)

    field = next(iter(conf.fields.values()))
    assert field.approx_order == 3


def test_override_material_dictionary_in_tuple():
    data = {'materials': {'m': ({'k': 1.0, 'density': 3.0},)}}
    override = {'materials': {'m': ({'k': 2.5},)}}
    saved = deepcopy(override)
    conf = ProblemConf.from_dict(data, None, required=[], override=override)

    material = next(iter(conf.materials.values()))
    assert material.values == {'k': 2.5, 'density': 3.0}
    assert override == saved


def test_override_boundary_tuple():
    data = {'ebcs': {'left': ('Left', {'u.all': 0.0})}}
    override = {'ebcs': {'left': ('Right', {'u.all': 3.25})}}
    conf = ProblemConf.from_dict(data, None, required=[], override=override)

    boundary = next(iter(conf.ebcs.values()))
    assert boundary.region == 'Right'
    assert boundary.dofs == {'u.all': 3.25}


def test_recursive_tuple_update_preserves_unrelated_items():
    data = {'solver': ('old', {'nested': {'keep': 1, 'replace': 2}}),
            'unrelated': 42}
    override = {'solver': ('new', {'nested': {'replace': 3, 'add': 4}})}
    saved = deepcopy(override)
    result = update_dict_recursively(data, override, tuples_too=True)

    assert result is data
    assert result == {
        'solver': ('new', {'nested': {'keep': 1, 'replace': 3, 'add': 4}}),
        'unrelated': 42,
    }
    assert override == saved
