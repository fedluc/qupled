import numpy as np

from dataclasses import dataclass
from typing import get_type_hints


def serializable_dataclass(cls):

    cls = dataclass(cls)

    def to_dict(self):
        result = {}
        for key, value in self.__dict__.items():
            if hasattr(value, "to_dict") and callable(value.to_dict):
                result[key] = value.to_dict()
            elif isinstance(value, np.ndarray):
                result[key] = value.tolist()
            else:
                result[key] = value
        return result

    @classmethod
    def from_dict(cls, d):
        obj = cls.__new__(cls)
        annotations = get_type_hints(cls)
        print("Deserializing object of type:", cls.__name__)
        print(annotations)
        for key, value in d.items():
            expected_type = annotations.get(key)
            from_dict_fn = getattr(expected_type, "from_dict", None)
            convert_to_np_array = expected_type is np.ndarray and isinstance(
                value, list
            )
            call_from_dict = callable(from_dict_fn) and isinstance(value, dict)
            if convert_to_np_array:
                setattr(obj, key, np.array(value))
            elif call_from_dict:
                setattr(obj, key, from_dict_fn(value))
            else:
                setattr(obj, key, value)
        return obj

    cls.to_dict = to_dict
    cls.from_dict = from_dict
    return cls
