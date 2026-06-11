#!/usr/bin/env python3

import os


def main():
    value = os.getenv("OPENSN_TEST_ENV_VALUE", "<unset>")
    number = os.getenv("OPENSN_TEST_ENV_NUMBER", "<unset>")
    print(f"ENV_VALUE {value}")
    print(f"ENV_NUMBER {number}")


if __name__ == "__main__":
    main()
