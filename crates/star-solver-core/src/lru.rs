//! 零依赖 LRU 缓存（HashMap + VecDeque 实现，容量固定）。

use std::collections::{HashMap, VecDeque};
use std::hash::Hash;

pub struct LruCache<K: Clone + Eq + Hash, V> {
    map: HashMap<K, V>,
    order: VecDeque<K>,
    cap: usize,
}

impl<K: Clone + Eq + Hash, V> LruCache<K, V> {
    pub fn new(cap: usize) -> Self {
        Self {
            map: HashMap::with_capacity(cap.max(1)),
            order: VecDeque::with_capacity(cap.max(1)),
            cap: cap.max(1),
        }
    }

    pub fn contains(&self, k: &K) -> bool {
        self.map.contains_key(k)
    }

    /// 读取并把键提到最新位置。
    pub fn get(&mut self, k: &K) -> Option<&V> {
        if self.map.contains_key(k) {
            if let Some(pos) = self.order.iter().position(|x| x == k) {
                let key = self.order.remove(pos).unwrap();
                self.order.push_back(key);
            }
            self.map.get(k)
        } else {
            None
        }
    }

    pub fn put(&mut self, k: K, v: V) {
        if self.map.contains_key(&k) {
            if let Some(pos) = self.order.iter().position(|x| x == &k) {
                let key = self.order.remove(pos).unwrap();
                self.order.push_back(key);
            }
            self.map.insert(k, v);
            return;
        }
        if self.order.len() >= self.cap {
            if let Some(old) = self.order.pop_front() {
                self.map.remove(&old);
            }
        }
        self.order.push_back(k.clone());
        self.map.insert(k, v);
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}
