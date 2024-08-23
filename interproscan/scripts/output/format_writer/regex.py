import re


NT_SEQ_ID_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=(\d+\.+\d+)\s+.+frame=(\d+)\s+desc=(.*)$")
NT_KEY_PATTERN = re.compile(r"^(.*)_orf\d+$")
