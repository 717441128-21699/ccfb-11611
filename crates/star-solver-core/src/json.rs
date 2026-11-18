//! 极简 JSON 解析器（仅满足 manifest 与结果输出所需的子集）。
//!
//! 支持 object / array / string / number / true / false / null，
//! 不带数字索引下标语法、不支持注释。核心库因此不需要 serde。

use std::collections::HashMap;

#[derive(Debug, Clone)]
pub enum JsonValue {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<JsonValue>),
    Obj(HashMap<String, JsonValue>),
}

impl JsonValue {
    pub fn as_object(&self) -> Option<&HashMap<String, JsonValue>> {
        match self {
            JsonValue::Obj(o) => Some(o),
            _ => None,
        }
    }
    pub fn as_array(&self) -> Option<&Vec<JsonValue>> {
        match self {
            JsonValue::Arr(a) => Some(a),
            _ => None,
        }
    }
    pub fn as_str(&self) -> Option<&str> {
        match self {
            JsonValue::Str(s) => Some(s),
            _ => None,
        }
    }
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            JsonValue::Num(n) => Some(*n),
            _ => None,
        }
    }
    pub fn as_u64(&self) -> Option<u64> {
        self.as_f64().map(|n| n as u64)
    }
}

/// object 访问缺键返回 Null，方便链式 `o["k"].as_f64()`。
impl std::ops::Index<&str> for JsonValue {
    type Output = JsonValue;
    fn index(&self, key: &str) -> &JsonValue {
        match self {
            JsonValue::Obj(o) => o.get(key).unwrap_or(&JsonValue::Null),
            _ => &JsonValue::Null,
        }
    }
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Parser<'a> {
    fn skip_ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn parse_value(&mut self) -> Result<JsonValue, String> {
        self.skip_ws();
        if self.i >= self.b.len() {
            return Err("unexpected eof".into());
        }
        Ok(match self.b[self.i] {
            b'{' => self.parse_object()?,
            b'[' => self.parse_array()?,
            b'"' => JsonValue::Str(self.parse_string()?),
            b't' | b'f' => self.parse_bool()?,
            b'n' => self.parse_null()?,
            b'-' | b'0'..=b'9' => self.parse_number()?,
            c => return Err(format!("unexpected char {c} at {}", self.i)),
        })
    }

    fn parse_object(&mut self) -> Result<JsonValue, String> {
        self.i += 1; // {
        let mut map = HashMap::new();
        self.skip_ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Ok(JsonValue::Obj(map));
        }
        loop {
            self.skip_ws();
            let key = self.parse_string()?;
            self.skip_ws();
            if self.next() != Some(b':') {
                return Err("expected ':'".into());
            }
            let val = self.parse_value()?;
            map.insert(key, val);
            self.skip_ws();
            match self.next() {
                Some(b',') => continue,
                Some(b'}') => break,
                other => return Err(format!("expected , or }} got {other:?}")),
            }
        }
        Ok(JsonValue::Obj(map))
    }

    fn parse_array(&mut self) -> Result<JsonValue, String> {
        self.i += 1; // [
        let mut arr = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Ok(JsonValue::Arr(arr));
        }
        loop {
            arr.push(self.parse_value()?);
            self.skip_ws();
            match self.next() {
                Some(b',') => continue,
                Some(b']') => break,
                other => return Err(format!("expected , or ] got {other:?}")),
            }
        }
        Ok(JsonValue::Arr(arr))
    }

    fn parse_string(&mut self) -> Result<String, String> {
        if self.next() != Some(b'"') {
            return Err("expected string".into());
        }
        let mut out = String::new();
        loop {
            match self.next() {
                None => return Err("unterminated string".into()),
                Some(b'"') => break,
                Some(b'\\') => match self.next() {
                    Some(b'"') => out.push('"'),
                    Some(b'\\') => out.push('\\'),
                    Some(b'/') => out.push('/'),
                    Some(b'n') => out.push('\n'),
                    Some(b't') => out.push('\t'),
                    Some(b'r') => out.push('\r'),
                    Some(b'b') => out.push('\u{0008}'),
                    Some(b'f') => out.push('\u{000C}'),
                    Some(b'u') => {
                        let mut code = 0u32;
                        for _ in 0..4 {
                            let c = self.next().ok_or("bad unicode escape")?;
                            code = code * 16
                                + (c as char).to_digit(16).ok_or("bad hex digit")?;
                        }
                        if let Some(ch) = char::from_u32(code) {
                            out.push(ch);
                        }
                    }
                    _ => return Err("bad escape".into()),
                },
                Some(c) => {
                    // 多字节 UTF-8 原样收集
                    if c < 0x80 {
                        out.push(c as char);
                    } else {
                        let start = self.i - 1;
                        let len = if c >= 0xF0 {
                            4
                        } else if c >= 0xE0 {
                            3
                        } else {
                            2
                        };
                        while self.i < start + len {
                            self.i += 1;
                        }
                        out.push_str(
                            std::str::from_utf8(&self.b[start..start + len])
                                .map_err(|_| "bad utf8")?,
                        );
                    }
                }
            }
        }
        Ok(out)
    }

    fn parse_bool(&mut self) -> Result<JsonValue, String> {
        if self.b[self.i..].starts_with(b"true") {
            self.i += 4;
            Ok(JsonValue::Bool(true))
        } else if self.b[self.i..].starts_with(b"false") {
            self.i += 5;
            Ok(JsonValue::Bool(false))
        } else {
            Err("bad literal".into())
        }
    }

    fn parse_null(&mut self) -> Result<JsonValue, String> {
        if self.b[self.i..].starts_with(b"null") {
            self.i += 4;
            Ok(JsonValue::Null)
        } else {
            Err("bad literal".into())
        }
    }

    fn parse_number(&mut self) -> Result<JsonValue, String> {
        let start = self.i;
        if self.peek() == Some(b'-') {
            self.i += 1;
        }
        while matches!(self.peek(), Some(c) if c.is_ascii_digit() || c == b'.' || c == b'e' || c == b'E' || c == b'+' || c == b'-')
        {
            self.i += 1;
        }
        std::str::from_utf8(&self.b[start..self.i])
            .map_err(|e| e.to_string())?
            .parse::<f64>()
            .map(JsonValue::Num)
            .map_err(|e| e.to_string())
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }
    fn next(&mut self) -> Option<u8> {
        let c = self.b.get(self.i).copied();
        if c.is_some() {
            self.i += 1;
        }
        c
    }
}

pub fn parse(text: &str) -> Result<JsonValue, String> {
    let mut p = Parser {
        b: text.as_bytes(),
        i: 0,
    };
    p.parse_value()
}

/// JSON 字符串转义。
pub fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}
