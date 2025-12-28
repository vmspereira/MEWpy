# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
File containing string expression parsing utility functions.
Expressions are decomposed into binary parsed trees.

Author: Vitor Pereira
##############################################################################
"""
import logging
import re
import sys
import typing as T
from abc import abstractmethod
from copy import copy
from math import (
    acos,
    asin,
    atan,
    atan2,
    ceil,
    cos,
    cosh,
    degrees,
    e,
    exp,
    floor,
    log,
    log2,
    log10,
    pi,
    radians,
    sin,
    sinh,
    sqrt,
    tan,
    tanh,
)
from operator import add, mul, pow, sub, truediv

logger = logging.getLogger(__name__)

# Boolean operator symbols
S_AND = "&"
S_OR = "|"
S_NOT = "~"
S_ON = "1"
S_OFF = "0"
S_GREATER = ">"
S_LESS = "<"
S_EQUAL = "="
S_GREATER_THAN_EQUAL = ">="
S_LESS_THAN_EQUAL = "<="

# Empty leaf symbol
EMPTY_LEAF = "@"

BOOLEAN_OPERATORS = [S_AND, S_OR, S_NOT]

BOOLEAN_STATES = [S_ON, S_OFF]

BOOLEAN_RELATIONALS = [S_GREATER, S_LESS, S_EQUAL]
BOOLEAN_EQUAL_RELATIONALS = [S_GREATER_THAN_EQUAL, S_LESS_THAN_EQUAL]

# Increases the system recursion limit for long expressions
sys.setrecursionlimit(100000)

# latex #############################################


class Latex:
    def __init__(self, text) -> None:
        """A simple class for _repr_latex_

        Args:
            text (str): LaTeX string
        """
        self.text = text

    def __repr__(self) -> str:
        return self.text

    def _repr_latex_(self) -> str:
        """Return LaTeX for Jupyter display.

        Jupyter needs $$ delimiters to recognize and render LaTeX as math.
        """
        return "$$ %s $$" % self.text


def _escape_latex_text(text: str) -> str:
    """Escape special LaTeX characters in text for use in \\text{} commands.

    Args:
        text: The text string to escape.

    Returns:
        The escaped text safe for LaTeX rendering.

    Note:
        This escapes characters that have special meaning in LaTeX:
        - Underscore (_) is used for subscripts
        - Caret (^) is used for superscripts
        - Braces ({}) are used for grouping
        - Percent (%) starts comments
        - Ampersand (&) is used in tables
        - Hash (#) is used for macro parameters
        - Dollar ($) enters/exits math mode
        - Tilde (~) is a non-breaking space
        - Backslash (\\) starts commands
    """
    # Order matters: escape backslash first, then other characters
    replacements = [
        ("\\", r"\textbackslash{}"),
        ("_", r"\_"),
        ("^", r"\^{}"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("%", r"\%"),
        ("&", r"\&"),
        ("#", r"\#"),
        ("$", r"\$"),
        ("~", r"\textasciitilde{}"),
    ]

    for char, replacement in replacements:
        text = text.replace(char, replacement)

    return text


def convert_constant(value: T.Any) -> str:
    """Helper to convert constant values to LaTeX string.

    Args:
        value: A constant value.

    Returns:
        The LaTeX representation of `value`.

    Note:
        String values are escaped to handle special LaTeX characters like
        underscores, which would otherwise cause rendering errors in notebooks.
    """
    if value is None or isinstance(value, bool):
        return r"\mathrm{" + str(value) + "}"
    if isinstance(value, (int, float, complex)):
        # TODO(odashi): Support other symbols for the imaginary unit than j.
        # Current implementation only supports 'j' for complex numbers
        # Consider using a configurable symbol or supporting both 'j' and 'i'
        return str(value)
    if isinstance(value, str):
        # Escape special LaTeX characters in strings (e.g., underscores in 'value_1')
        escaped_value = _escape_latex_text(value)
        return r"\textrm{" + escaped_value + "}"
    if isinstance(value, bytes):
        escaped_value = _escape_latex_text(str(value))
        return r"\textrm{" + escaped_value + "}"
    if value is ...:
        return r"\cdots"
    raise ValueError(f"Unrecognized constant: {type(value).__name__}")


def paren(src: str) -> str:
    """Adds surrounding parentheses: "(" and ")"."""
    return r"\mathopen{}\left( " + src + r" \mathclose{}\right)"


# TODO: Add more functions
latex = {
    "*": lambda x, y: rf"{x} \times {y}",
    "/": lambda x, y: r"\frac {" + x + "} {" + y + "}",
    "+": lambda x, y: f"{x} + {y}",
    "-": lambda x, y: f"{x} - {y}",
    "pow": lambda x, y: r"{" + x + r"}^{" + y + r"}",
    "^": lambda x, y: r"{" + x + r"}^{" + y + r"}",
    "sqrt": lambda x, y: r"\sqrt {" + y + r"}",
}

# Operators precedence used to add parentheses when
# needed as they are removed in the parsing tree
MAX_PRECEDENCE = 10
latex_precedence = {
    "+": 0,
    "-": 0,
    "*": 1,
    "/": 1,
    "^": 2,
    "pow": 2,
}


# Evaluate #############################################


def evaluate_expression(expression: str, variables: T.List[str]) -> T.Any:
    """Evaluates a logical expression (containing variables, 'and','or','(' and ')')
    against the presence (True) or absence (False) of propositions within a list.
    The evaluation is achieved using a safe parsing tree approach.

    :param str expression: The expression to be evaluated.
    :param list variables: List of variables to be evaluated as True.
    :returns: A boolean evaluation of the expression.

    """
    # Use the tree-based evaluator which is safer than eval()
    return evaluate_expression_tree(expression, variables)


def evaluate_expression_tree(expression: str, variables: T.List[str]) -> T.Any:
    """Evaluates a logical expression (containing variables,
    'and','or', 'not','(' , ')') against the presence (True)
    or absence (False) of propositions within a list.
    Assumes the correctness of the expression. The evaluation
    is achieved by means of a parsing tree.

    :param str expression: The expression to be evaluated.
    :param list variables: List of variables to be evaluated as True.
    :returns: A boolean evaluation of the expression.

    """
    t = build_tree(expression, Boolean)
    evaluator = BooleanEvaluator(variables)
    res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
    return res


def maybe_fn(f: T.Callable, v1: T.Any, v2: T.Any) -> T.Any:
    """Maybe evaluator: if one of the arguments is None, it
    retuns the value of the other argument. If both arguments
    are None, it returns None. If both arguments are not None
    it returns the evaluation f(v1,v2).

    :param f: a function
    :param v1: the first argument
    :param v2: the second argument
    """
    if v1 is None:
        return v2
    elif v2 is None:
        return v1
    else:
        return f(v1, v2)


# Parsing Tree #############################################


class Node(object):
    """
    Binary syntax tree node.

    :param value: The node value.
    :param left: The left node or None.
    :param right: The right node or None.
    """

    def __init__(
        self,
        value: T.Any,
        left: T.Union["Node", None] = None,
        right: T.Union["Node", None] = None,
        tp=0,
    ) -> None:
        """Binary tree.  Empty leafs with None or 'EMPTY_LEAF'
           are always left sided.

        Args:
            value : Node value
            left (Node, optional): [description]. Defaults to None.
            right (Node, optional): [description]. Defaults to None.
            tp (int, optional): 1-Unary function, 2-binary function. Defaults to 0.

        Raises:
            ValueError
        """
        if left is not None and not isinstance(left, Node):
            raise ValueError("Invalid right element")
        if right is not None and not isinstance(right, Node):
            raise ValueError("Invalid right element")
        self.value = value
        self.left = left
        self.right = right
        self.tp = tp

    def __repr__(self) -> str:
        if self.is_leaf():
            return str(self.value)
        else:
            return f"{str(self.value)} " f"( {str(self.left)} ," f" {str(self.right)} )"

    # def _repr_latex_(self):
    #    return "$$ %s $$" % (self.to_latex())

    def __str__(self) -> str:
        return self.to_infix()

    def is_leaf(self) -> bool:
        """
        :returns: True if the node is a leaf False otherwise.
        Both left and right are None.
        """
        return not self.left and not self.right

    def is_empty_leaf(self) -> bool:
        """
        :returns: True if the node is an empty leaf False otherwise.
        """
        return self.value == EMPTY_LEAF

    def is_unary(self) -> bool:
        """
        :returns: True if the node is a unary operation False otherwise.
        """
        return (self.left.is_empty_leaf() and not self.right.is_empty_leaf()) or (
            not self.left.is_empty_leaf() and self.right.is_empty_leaf()
        )

    def is_binary(self) -> bool:
        """
        :returns: True if the node is a binary operation False otherwise.
        """
        return not self.left.is_empty_leaf() and not self.right.is_empty_leaf()

    def get_operands(self):
        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return set()
            else:
                return {self.value}
        else:
            return self.left.get_operands().union(self.right.get_operands())

    def get_operators(self):
        if self.is_leaf():
            return set()
        else:
            return {self.value}.union(self.left.get_operators()).union(self.right.get_operators())

    def get_parameters(self):
        """Parameters are all non numeric symbols in an expression"""
        if self.is_leaf():
            if self.value != EMPTY_LEAF and not is_number(self.value):
                return {self.value}
            else:
                return set()
        else:
            return self.left.get_parameters().union(self.right.get_parameters())

    def print_node(self, level=0):
        """Prints a parsing tree of the expression."""
        tabs = ""
        for _ in range(level):
            tabs += "\t"
        if self.is_leaf():
            if self.value != EMPTY_LEAF:
                logger.debug(f"{tabs}|____{self.value}")
            else:
                pass
        else:
            logger.debug(f"{tabs}|____{self.value}")
        if self.left is not None:
            self.left.print_node(level + 1)
        if self.right is not None:
            self.right.print_node(level + 1)

    def evaluate(self, f_operand=None, f_operator=None):
        """
        Evaluates the expression using the f_operand and
        f_operator mapping functions
        """
        if f_operand is None or f_operator is None:
            raise ValueError(
                "Both f_operand and f_operator functions must be provided for safe evaluation. "
                "Using eval() has been disabled for security reasons. "
                "Please provide appropriate evaluator functions."
            )
        elif self.is_leaf():
            return f_operand(self.value)
        else:
            return maybe_fn(
                f_operator(self.value),
                self.left.evaluate(f_operand, f_operator),
                self.right.evaluate(f_operand, f_operator),
            )

    def get_conditions(self):
        """
        Retrieves the propositional conditions

        return: The list of conditions
        """
        ops = self.get_operands()
        return set([i for i in ops if is_condition(i)])

    def replace(self, r_map: dict):
        """Apply the mapping replacing to the tree

        Args:
            map (dict): replacement mapping

        Returns:
            Node: new Tree with replace entries
        """
        if self.is_leaf():
            v = r_map[self.value] if self.value in r_map.keys() else self.value
            return Node(v, None, None)
        else:
            return Node(self.value, self.left.replace(r_map), self.right.replace(r_map), self.tp)

    def replace_node(self, value, node):

        if self.value is not None and self.value == value:
            self.value = node.value
            self.left = node.left.copy()
            self.right = node.right.copy()
            self.tp = node.tp

        elif not self.is_leaf():
            self.left.replace_node(value, node)
            self.right.replace_node(value, node)

        else:
            pass

    def replace_nodes(self, nodes: dict):
        for k, v in nodes.items():
            self.replace_node(k, v)

    def to_infix(
        self,
        opar: str = "(",
        cpar: str = ")",
        sep: str = " ",
        fsep: str = " , ",
        replacers=None,
    ) -> str:
        """Infix string representation

        :param opar: open parentheses string, defaults to '( '
        :type opar: str, optional
        :param cpar: close parentheses string, defaults to ' )'
        :type cpar: str, optional
        :param sep: symbols separator, defaults to ' '
        :type sep: str, optional
        :param fsep: function argument separator, defaults to ' , '
        :type fsep: str, optional
        :return: An infix string representation of the node
        :rtype: str
        """

        rep = {S_AND: "and", S_OR: "or", "^": "**"}
        if replacers:
            rep.update(replacers)

        def rval(value):
            return str(rep[value]) if value in rep.keys() else str(value)

        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return ""
            else:
                return rval(self.value)
        elif self.tp >= 2:
            op = opar if self.tp == 2 else ""
            cp = cpar if self.tp == 2 else ""
            return "".join(
                [
                    rval(self.value),
                    opar,
                    self.left.to_infix(op, cp, sep, fsep),
                    fsep,
                    self.right.to_infix(op, cp, sep, fsep),
                    cpar,
                ]
            )
        elif self.tp == 1:
            return "".join(
                [
                    rval(self.value),
                    opar,
                    self.right.to_infix(opar, cpar, sep, fsep),
                    cpar,
                ]
            )
        else:
            return "".join(
                [
                    opar,
                    self.left.to_infix(opar, cpar, sep, fsep),
                    sep,
                    rval(self.value),
                    sep,
                    self.right.to_infix(opar, cpar, sep, fsep),
                    cpar,
                ]
            )

    def to_latex(self) -> T.Tuple[str, int]:
        """
        Simple conversion of a parsing tree to LaTeX.
        Operator precedences are used to decide when to add
        parentheses.

        Returns:
            LaTeX str, precedence of last operator
        """
        if self.is_leaf():
            if self.value == EMPTY_LEAF:
                return "", MAX_PRECEDENCE
            else:
                return convert_constant(self.value), MAX_PRECEDENCE
        else:
            op = self.value.strip()
            if op in latex:
                left_latex, left_prec = self.left.to_latex()
                right_latex, right_prec = self.right.to_latex()
                p = latex_precedence.get(op, MAX_PRECEDENCE)
                s_l = paren(left_latex) if p > left_prec else left_latex
                s_r = paren(right_latex) if p > right_prec else right_latex
                return latex[op](s_l, s_r), p

            elif self.tp == 1:
                return (
                    convert_constant(self.value) + r"\left(" + self.right.to_latex()[0] + r"\right)",
                    MAX_PRECEDENCE,
                )
            else:
                return (
                    convert_constant(self.value)
                    + r"\left("
                    + self.left.to_latex()[0]
                    + ","
                    + self.right.to_latex()[0]
                    + r"\right)",
                    MAX_PRECEDENCE,
                )

    def copy(self):
        if self.is_leaf():
            return Node(copy(self.value), None, None)
        else:
            return Node(copy(self.value), self.left.copy(), self.right.copy(), self.tp)


class Syntax:
    """Defines an interface for the tree syntax parsing
    with operators, their precedence and associativity.
    """

    operators = []

    @staticmethod
    @abstractmethod
    def is_operator(op):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def is_greater_precedence(op1, op2):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def associativity(op):
        raise NotImplementedError

    @staticmethod
    def arity(op):
        return 2

    @staticmethod
    def replace():
        return {}

    @staticmethod
    def sub(op):
        return op


class Arithmetic(Syntax):
    """Defines a basic arithmetic syntax."""

    operators = ["+", "-", "**", "*", "/", "^"]

    @staticmethod
    def is_operator(op):
        return op in ["+", "-", "*", "/", "^"]

    @staticmethod
    def is_greater_precedence(op1, op2):
        pre = {"+": 0, "-": 0, "*": 1, "/": 1, "^": 2}
        return pre[op1] >= pre[op2]

    @staticmethod
    def associativity(op):
        ass = {"+": 0, "-": 0, "*": 0, "/": 0, "^": 1}
        return ass[op]

    @staticmethod
    def arity(op):
        ar = {"+": 2, "-": 2, "*": 2, "/": 2, "^": 2}
        return ar[op]

    @staticmethod
    def sub(op):
        if op == "**":
            return "^"
        else:
            return op

    @staticmethod
    def rsub(op):
        if op == "^":
            return "**"
        else:
            return op


class ArithmeticEvaluator:
    @staticmethod
    def f_operator(op):
        operators = {"+": add, "-": sub, "*": mul, "/": truediv, "^": pow, "pow": pow}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    @staticmethod
    def f_operand(op):
        if is_number(op):
            return float(op)
        else:
            raise ValueError(f"Operator {op} not defined")


class Boolean(Syntax):
    """
    A boolean syntax parser where (NOT) is considered to
    be defined as a binary operator (NOT a) = (EMPTY NOT a)
    """

    operators = [S_AND, S_OR, S_NOT]

    @staticmethod
    def is_operator(op):
        return op in [S_AND, S_OR, S_NOT]

    @staticmethod
    def is_greater_precedence(op1, op2):
        pre = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return pre[op1] >= pre[op2]

    @staticmethod
    def associativity(op):
        ass = {S_AND: 0, S_OR: 0, S_NOT: 1}
        return ass[op]

    @staticmethod
    def arity(op):
        ar = {S_AND: 2, S_OR: 2, S_NOT: 1}
        return ar[op]

    @staticmethod
    def replace():
        r = {
            "not": [EMPTY_LEAF, S_NOT],
            "and": [S_AND],
            "or": [S_OR],
            "&": [S_AND],
            "|": [S_OR],
            "~": [EMPTY_LEAF, S_NOT],
        }
        return r


class BooleanEvaluator:
    """A boolean evaluator.

    :param list true_list: Operands evaluated as True. Operands not in the list are evaluated as False
    :param dict variables: A dictionary mapping symbols to values. Used to evaluate conditions.

    """

    def __init__(self, true_list=None, variables=None):
        self.true_list = true_list if true_list is not None else []
        self.vars = variables if variables is not None else {}

    def f_operator(self, op):
        operators = {
            S_AND: lambda x, y: x and y,
            S_OR: lambda x, y: x or y,
            S_NOT: lambda x, y: not y,
        }
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self, op):
        if op.upper() == "TRUE" or op == "1" or op in self.true_list:
            return True
        elif is_condition(op):
            # Use safe condition evaluator instead of eval()
            try:
                return evaluate_condition(op, self.vars)
            except ValueError:
                # If condition evaluation fails (invalid format or non-numeric values), default to False
                # This can happen with malformed conditions like "x y z" or "x > abc"
                return False
        else:
            return False

    def set_true_list(self, true_list):
        self.true_list = true_list


class GeneEvaluator:
    """An evaluator for genes expression.

    :param genes_value: A dictionary mapping genes to values. Gene not in the dictionary have a value of 1.
    :param and_operator: function to be applied instead of (and', '&') (not case sensitive)
    :param or_operator: function to be applied instead of ('or','|') (not case sensitive)
    """

    def __init__(
        self,
        genes_value,
        and_operator=min,
        or_operator=max,
        prefix="",
        unexpressed_value=1,
    ):

        self.genes_value = genes_value
        self.and_operator = and_operator
        self.or_operator = or_operator
        self.prefix = prefix
        self.unexpressed_value = unexpressed_value

    def f_operator(self, op):
        operators = {S_AND: self.and_operator, S_OR: self.or_operator}
        if op in operators.keys():
            return operators[op]
        else:
            raise ValueError(f"Operator {op} not defined")

    def f_operand(self, op):
        if op[len(self.prefix) :] in self.genes_value:
            return self.genes_value[op]
        else:
            return self.unexpressed_value


def tokenize_function(exp: str) -> T.List[str]:
    """Tokenize a function "f(...)" expression string

    Args:
        exp (str): expression

    Returns:
        T.List[str]: List of tokens
    """
    p = 0
    s = -1
    tokens = []
    i = 0
    while i < len(exp):
        if exp[i] == "(":
            s += 1
            if s == 0:
                p = i + 1
            if not tokens:
                tokens = [exp[:i].strip()]
        elif exp[i] == ")":
            s -= 1
            if s == -1:
                tokens.append(exp[p:i].strip())
        elif exp[i] == "," and s == 0:
            tokens.append(exp[p:i].strip())
            p = i + 1
        i += 1
    if not tokens:
        return [exp]
    else:
        return tokens


def list2tree(values, rules):
    if len(values) == 0:
        return Node(EMPTY_LEAF)
    elif len(values) == 1:
        return build_tree(values[0], rules)
    else:
        return Node(",", build_tree(values[0], rules), list2tree(values[1:], rules))


# Tree
def build_tree(exp: str, rules: Syntax) -> Node:
    """
    Builds a parsing syntax tree for basic mathematical expressions

    :param exp: the expression to be parsed
    :param rules: Syntax definition rules
    """
    assert exp.count("(") == exp.count(")"), "The expression is parentheses unbalanced."
    replace_dic = rules.replace()
    exp_ = tokenize_infix_expression(exp, rules)
    exp_list = []
    for token in exp_:
        if token.lower() in replace_dic:
            exp_list.extend(replace_dic[token.lower()])
        else:
            exp_list.append(token)
    stack = []
    tree_stack = []
    predecessor = None
    i = 0
    while i < len(exp_list):
        token = exp_list[i]
        if not (rules.is_operator(token) or token in ["(", ")"]):
            if i < len(exp_list) - 2 and exp_list[i + 1] == "(":
                s = 1
                p = i + 2
                while p < len(exp_list) and s > 0:
                    if exp_list[p] == "(":
                        s += 1
                    elif exp_list[p] == ")":
                        s -= 1
                    p += 1
                token = " ".join(exp_list[i:p])
                i = p - 1

            if predecessor and not (rules.is_operator(predecessor) or predecessor in ["(", ")"]):
                s = tree_stack[-1].value
                tree_stack[-1].value = s + " " + token
            else:
                if "(" in token:
                    f = tokenize_function(token)
                    fname = f[0]
                    params = f[1:]
                    if len(params) == 1:
                        t = Node(fname, Node(EMPTY_LEAF), build_tree(params[0], rules), 1)
                    elif len(params) == 2:
                        t = Node(fname, build_tree(params[0], rules), build_tree(params[1], rules), 2)
                    else:
                        t = Node(fname, build_tree(params[0], rules), list2tree(params[1:], rules), len(f) - 1)
                else:
                    t = Node(token)
                tree_stack.append(t)
        elif rules.is_operator(token):
            if not stack or stack[-1] == "(":
                stack.append(token)

            elif rules.is_greater_precedence(token, stack[-1]) and rules.associativity(token) == 1:
                stack.append(token)

            else:
                while (
                    stack
                    and stack[-1] != "("
                    and rules.is_greater_precedence(stack[-1], token)
                    and rules.associativity(token) == 0
                ):
                    popped_item = stack.pop()
                    t = Node(popped_item)
                    t1 = tree_stack.pop()
                    t2 = tree_stack.pop()
                    t.right = t1
                    t.left = t2
                    tree_stack.append(t)
                stack.append(token)

        elif token == "(":
            stack.append("(")

        elif token == ")":
            while stack[-1] != "(":
                popped_item = stack.pop()
                t = Node(popped_item)
                t1 = tree_stack.pop()
                t2 = tree_stack.pop()
                t.right = t1
                t.left = t2
                tree_stack.append(t)
            stack.pop()
        predecessor = token
        i += 1

    while stack:
        popped_item = stack.pop()
        t = Node(popped_item)
        t1 = tree_stack.pop()
        t2 = tree_stack.pop()
        t.right = t1
        t.left = t2
        tree_stack.append(t)

    t = tree_stack.pop()

    return t


def tokenize_infix_expression(exp: str, rules: Syntax = None) -> T.List[str]:
    _exp = exp.replace("(", " ( ").replace(")", " ) ")
    if rules:
        for op in rules.operators:
            _exp = _exp.replace(op, " " + rules.sub(op) + " ")
    tokens = _exp.split(" ")
    return list(filter(lambda x: x != "", tokens))


def is_number(token: str) -> bool:
    """Returns True if the token is a number"""
    return token.replace(".", "", 1).replace("-", "", 1).isnumeric()


def is_condition(token: str) -> bool:
    """Returns True if the token is a condition"""
    regexp = re.compile(r">|<|=")
    return bool(regexp.search(token))


def evaluate_condition(condition: str, variables: T.Dict[str, T.Union[int, float]]) -> bool:
    """Safely evaluates a condition expression like 'x > 5' or 'y <= 10'.

    This function parses and evaluates comparison expressions without using eval(),
    providing a safer alternative for regulatory condition evaluation.

    :param condition: The condition string (e.g., 'x > 5', 'y <= 10')
    :param variables: Dictionary mapping variable names to numeric values
    :returns: Boolean result of the condition evaluation
    :raises ValueError: If the condition format is invalid
    """
    condition = condition.strip()

    # Define comparison operators
    operators = {
        ">=": lambda x, y: x >= y,
        "<=": lambda x, y: x <= y,
        "=>": lambda x, y: x >= y,  # Alternative notation
        "=<": lambda x, y: x <= y,  # Alternative notation
        ">": lambda x, y: x > y,
        "<": lambda x, y: x < y,
        "==": lambda x, y: x == y,
        "=": lambda x, y: x == y,
        "!=": lambda x, y: x != y,
    }

    # Try to parse the condition with different operator patterns
    # Sort operators by length (longest first) to match '>=' before '>'
    for op_str in sorted(operators.keys(), key=len, reverse=True):
        if op_str in condition:
            parts = condition.split(op_str, 1)
            if len(parts) == 2:
                left_str, right_str = parts[0].strip(), parts[1].strip()

                # Determine which side is the variable and which is the value
                left_is_var = not is_number(left_str)
                right_is_var = not is_number(right_str)

                if left_is_var and not right_is_var:
                    # Format: var op value (e.g., 'x > 5')
                    var_name = left_str
                    value = float(right_str)
                    var_value = variables.get(var_name, 0)
                    return operators[op_str](var_value, value)

                elif not left_is_var and right_is_var:
                    # Format: value op var (e.g., '5 < x')
                    value = float(left_str)
                    var_name = right_str
                    var_value = variables.get(var_name, 0)
                    return operators[op_str](value, var_value)

                elif left_is_var and right_is_var:
                    # Format: var op var (e.g., 'x > y')
                    left_val = variables.get(left_str, 0)
                    right_val = variables.get(right_str, 0)
                    return operators[op_str](left_val, right_val)
                else:
                    # Format: value op value (e.g., '5 > 3')
                    left_val = float(left_str)
                    right_val = float(right_str)
                    return operators[op_str](left_val, right_val)

    # If no operator found, raise an error
    raise ValueError(f"Invalid condition format: '{condition}'")


def isozymes(exp: str) -> T.List[str]:
    """
    Parses a GPR and splits it into its isozymes as a list of strings.
    """
    tree = build_tree(exp, Boolean())

    def split_or(node):
        if node.is_leaf():
            return [node]
        elif node.is_binary():
            if node.value == S_AND:
                return [node]
            elif node.value == S_OR:
                return split_or(node.left) + split_or(node.right)
            else:
                raise ValueError(f"Unrecognized operator for node {node}")
        else:
            raise ValueError(f"{node} is not binary of leaf")

    prots = split_or(tree)

    # validate
    if not all([len(node.get_operators() - set(S_AND)) == 0 for node in prots]):
        raise ValueError(f"{exp} is a malformed expression")

    proteins = [node.to_infix(opar="", cpar="") for node in prots]
    return proteins
