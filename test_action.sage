from sage.categories.action import Action
from sage.structure.element import ModuleElement
from sage.modules.module import Module
import operator

class TwoVectors(Module):
    def __init__(self):
        Module.__init__(self, base = QQ)
        self.Element = TwoVectorsElement
        self.register_action(MatAction(MatrixSpace(QQ,2,2),self))
        self._populate_coercion_lists_()

    def _element_constructor_(self, x):
        return TwoVectorsElement(self, x)

class TwoVectorsElement(ModuleElement):
    def __init__(self, parent, x):
        ModuleElement.__init__(self, parent)
        self.val = (QQ**2)(x)
    def __repr__(self):
        return str(self.val)

class MatAction(Action):
    def __init__(self,G,M):
        Action.__init__(self,G,M,is_left = True,op = operator.mul)

    def _call_(self,g,v):
        return v.parent()(g * v.val)

V = TwoVectors()
v = V((1,2))
g = Matrix(QQ,2,2,[1,2,3,4])
w = g * v
